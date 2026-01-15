// ============================================================================
// MPM + Phase-Field (Allen–Cahn) + Heat(latent) + Gas convection
// - Uses quadgrid_cpp.h + particles.h
// - Reads input via DATA data("DATA.json") like your code #2 (mpm_data.h)
// - Mechanics: hyperelastic MPM with F/J/P + soften(phi)
// - Thermo: Allen–Cahn phase-field coupled to T with latent source rho*L*dphi/dt
// ============================================================================

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "quadgrid_cpp.h"
#include "particles.h"
#include "mpm_data.h"   // <-- uses DATA like code #2

using namespace std;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

struct Vec3 { double x=0,y=0,z=0; };

static inline int clampi(int a,int lo,int hi){ return max(lo,min(a,hi)); }
static inline double clampd(double x,double lo,double hi){ return max(lo,min(x,hi)); }
static inline double det2(double a00,double a01,double a10,double a11){ return a00*a11 - a01*a10; }
static inline int I2(int i,int j,int NX){ return j*NX+i; }

static inline void invT2(double F00,double F01,double F10,double F11,
                         double &iT00,double &iT01,double &iT10,double &iT11){
  double J=det2(F00,F01,F10,F11);
  double invJ=1.0/max(1e-12,J);
  iT00= F11*invJ;  iT01=-F01*invJ;
  iT10=-F10*invJ;  iT11= F00*invJ;
}

// Phase-field helpers
static inline double Wprime(double phi){
  // W = phi^2(1-phi)^2  => W' = 2phi(1-phi)(1-2phi)
  return 2.0*phi*(1.0-phi)*(1.0-2.0*phi);
}
static inline double gprime(double phi){
  // g = phi^2(3-2phi) => g' = 6phi(1-phi)
  return 6.0*phi*(1.0-phi);
}

struct Params {
  // Grid
  int NX=192, NY=192;
  double dx=1e-5, dy=1e-5;
  double thickness=1e-5;

  // Time control
  double dt=1e-10;
  int steps=200000;
  int output_every=200;
  int rebuild_every=10;

  // Gas thermal
  double rho_g=0.2;
  double cp_g=1500.0;
  double k_g=1.0;
  double T_inf=15000.0;
  Vec3 u_gas{300.0,0.0,0.0};
  bool advect=true;
  bool inflow_hot=true;
  bool dirichlet_gas=false;

  // Solid (Cr2O3)
  double rho_cr=5200.0;
  double cp_s=800.0;
  double cp_l=1000.0;
  double k_s=5.0;
  double k_l=3.0;
  double Tm=2435.0;
  double Lm=4e5; // J/kg

  // Mechanics
  double mu_s=1e8;
  double bulk_ratio=50.0;
  double soften_p=3.0;
  double mu_min=1e4;
  double K_min=1e5;
  double Jmin=0.60;
  double Jmax=1.40;
  double Lmax=2e6;

  // PIC blend & damping
  double pic_alpha=0.2;
  bool node_damp=true;
  double damp=0.99;

  // Initial shear
  bool use_initial_shear=true;
  double shear_v0=2.0;

  // Convection coupling (gas<->solid)
  bool use_conv=true;
  double qconv_max_Wm3=1e9;

  // h(U_rel)
  double h_ref=100.0;
  double U_ref=100.0;
  double h_exp=0.8;
  double h_min=10.0;
  double h_max=5e3;
  double Urel_max=1000.0;

  // Phase-field (Allen–Cahn)
  double M_phi      = 2e-6;
  double eps_phi    = 2.0e-6;
  double lambda_phi = 5e-4;
  double phi_min=0.0;
  double phi_max=1.0;

  string out_prefix="mpm2d_pf";
};

struct Cell{
  double T=300.0;
  double rhoCp=0.0;
  double k=0.0;
  double phi=0.0;   // phase-field order parameter
  double mask=0.0;  // occupancy [0..1]
};

struct Sim {
  using QGrid = quadgrid_t<std::vector<double>>;
  using Part  = particles_t;

  Params P;
  vector<Cell> grid;
  QGrid qg;
  unique_ptr<Part> mpm;

  map<string, vector<double>> vars;
  int nn=0;

  Sim(const Params& p):P(p){
    grid.resize((size_t)P.NX*P.NY);
    qg.set_sizes((QGrid::idx_t)P.NY,(QGrid::idx_t)P.NX,P.dy,P.dx);
    nn=(int)qg.num_local_nodes();

    auto make=[&](const string& k){ vars[k]=vector<double>(nn,0.0); };
    make("m");
    make("mom_vx"); make("mom_vy");
    make("F_ext_vx"); make("F_ext_vy");
    make("F_int_vx"); make("F_int_vy");
    make("Ftot_vx"); make("Ftot_vy");
    make("avx"); make("avy");
    make("vvx"); make("vvy");
  }

  static string pad5(int s){
    ostringstream oss; oss<<setw(5)<<setfill('0')<<s; return oss.str();
  }

  void clear_nodal(){
    for(auto& kv: vars) fill(kv.second.begin(), kv.second.end(), 0.0);
  }

  Vec3 cell_velocity(const Cell& c) const {
    if(c.mask>0.0) return Vec3{0,0,0};
    return P.u_gas;
  }

  void update_cell_properties(Cell& c){
    if(c.mask>0.0){
      double f=c.mask;
      // use phi as proxy to interpolate solid/liquid properties
      double rhoCp_cr = P.rho_cr*((1.0-c.phi)*P.cp_s + c.phi*P.cp_l);
      double k_cr     = (1.0-c.phi)*P.k_s + c.phi*P.k_l;
      c.rhoCp = f*rhoCp_cr + (1.0-f)*(P.rho_g*P.cp_g);
      c.k     = f*k_cr     + (1.0-f)*P.k_g;
    }else{
      c.rhoCp = P.rho_g*P.cp_g;
      c.k     = P.k_g;
      c.phi   = 0.0;
    }
  }

  // --------------------------------------------------------------------------
  // Particle seeding from DATA (like code2) or fallback disk if DATA lacks arrays
  // --------------------------------------------------------------------------
  void seed_particles_disk_from_DATA(const DATA& data){
  // Geometria disco da JSON
  Vec3 C{data.center.x, data.center.y, data.center.z};
  double R = data.R;
  int nP = data.nParticles;

  mt19937 rng(42);
  uniform_real_distribution<double> U(0.0,1.0);

  vector<double> xv(nP), yv(nP);

  double Adisk = M_PI*R*R;
  double Apt   = Adisk / (double)nP;
  double Vpt   = Apt * P.thickness;
  double mpt   = P.rho_cr * Vpt;

  for(int n=0;n<nP;++n){
    double r = R*sqrt(U(rng));
    double a = 2.0*M_PI*U(rng);
    xv[n] = C.x + r*cos(a);
    yv[n] = C.y + r*sin(a);
  }

  vector<string> iprops={};
  vector<string> dprops={
    "mp","V0",
    "vp_x","vp_y",
    "px","py",
    "F00","F01","F10","F11","J",
    "P00","P01","P10","P11",
    "dvx_dx","dvx_dy","dvy_dx","dvy_dy",
    "T","phi",
    "vp_x_old","vp_y_old",
    "Vp_scale"
  };

  mpm = make_unique<Part>((Part::idx_t)nP, iprops, dprops, qg, xv, yv);
  mpm->init_particle_mesh();
  mpm->build_mass();

  for(Part::idx_t p=0;p<mpm->num_particles;++p){
    mpm->dp("mp",p)=mpt;
    mpm->dp("V0",p)=Vpt;

    mpm->dp("F00",p)=1.0; mpm->dp("F01",p)=0.0;
    mpm->dp("F10",p)=0.0; mpm->dp("F11",p)=1.0;
    mpm->dp("J",p)=1.0;

    mpm->dp("T",p)=300.0;
    mpm->dp("phi",p)=0.0;

    if(P.use_initial_shear){
      double yy=mpm->y[p];
      // shear attorno al centro del disco
      mpm->dp("vp_x",p)= P.shear_v0*((yy-C.y)/max(1e-12,R));
      mpm->dp("vp_y",p)= 0.0;
    }else{
      mpm->dp("vp_x",p)=0.0;
      mpm->dp("vp_y",p)=0.0;
    }
  }
}


  // --------------------------------------------------------------------------
  // Rebuild mask + (optionally) T,phi from particles (like code1)
  // --------------------------------------------------------------------------
  void rebuild_mask_phi_T_from_particles(){
    const double Vcell=P.dx*P.dy*P.thickness;

    vector<double> mcell((size_t)P.NX*P.NY,0.0);
    vector<double> mp_sum((size_t)P.NX*P.NY,0.0);
    vector<double> mpT_sum((size_t)P.NX*P.NY,0.0);
    vector<double> mpPhi_sum((size_t)P.NX*P.NY,0.0);

    for(auto& c: grid){ c.mask=0.0; }

    for(Part::idx_t p=0;p<mpm->num_particles;++p){
      int i=clampi((int)floor(mpm->x[p]/P.dx),0,P.NX-1);
      int j=clampi((int)floor(mpm->y[p]/P.dy),0,P.NY-1);
      int id=I2(i,j,P.NX);

      double mp=mpm->dp("mp",p);
      mcell[id]     += mp;
      mp_sum[id]    += mp;
      mpT_sum[id]   += mp * mpm->dp("T",p);
      mpPhi_sum[id] += mp * mpm->dp("phi",p);
    }

    for(int j=0;j<P.NY;++j){
      for(int i=0;i<P.NX;++i){
        int id=I2(i,j,P.NX);
        Cell& c=grid[id];

        double f = mcell[id]/max(1e-30, P.rho_cr*Vcell);
        c.mask = clampd(f,0.0,1.0);

        if(c.mask>0.0 && mp_sum[id]>0.0){
          // optional: keep these as "sampling" from particles
          c.T   = mpT_sum[id]/mp_sum[id];
          c.phi = clampd(mpPhi_sum[id]/mp_sum[id],0.0,1.0);
        }else{
          c.phi=0.0;
          if(P.dirichlet_gas) c.T=300.0;
        }
        update_cell_properties(c);
      }
    }
  }

  void init_thermal_fields_once(){
    for(auto& c: grid){
      c.mask=0.0;
      c.phi=0.0;
      c.T=300.0;
      c.rhoCp=P.rho_g*P.cp_g;
      c.k=P.k_g;
    }
    rebuild_mask_phi_T_from_particles();
  }

  // --------------------------------------------------------------------------
  // h_eff from relative velocity (like code1)
  // --------------------------------------------------------------------------
  double compute_h_eff_from_relative_velocity(double &Urel_out) const {
    double msum=0.0, vxm=0.0;
    for(Part::idx_t p=0;p<mpm->num_particles;++p){
      double mp=mpm->dp("mp",p);
      double vx=clampd(mpm->dp("vp_x",p), -2000.0, 2000.0);
      msum += mp;
      vxm  += mp*vx;
    }
    double vx_s=(msum>0.0)?(vxm/msum):0.0;
    double urx = P.u_gas.x - vx_s;
    double Urel=fabs(urx);
    Urel=clampd(Urel,0.0,P.Urel_max);
    Urel_out=Urel;

    double ratio=Urel/max(1e-12,P.U_ref);
    double h=P.h_ref*pow(max(0.0,ratio),P.h_exp);
    return clampd(h,P.h_min,P.h_max);
  }

  // --------------------------------------------------------------------------
  // MPM mechanics (step style like code2)
  // --------------------------------------------------------------------------
  void step0_connectivity(){
    mpm->init_particle_mesh();
  }

  void step1_p2g_mass_momentum(){
    mpm->p2g(vars, {"mp"}, {"m"}, false);

    for(Part::idx_t p=0;p<mpm->num_particles;++p){
      mpm->dp("px",p)=mpm->dp("mp",p)*mpm->dp("vp_x",p);
      mpm->dp("py",p)=mpm->dp("mp",p)*mpm->dp("vp_y",p);
    }
    mpm->p2g(vars, {"px","py"}, {"mom_vx","mom_vy"}, false);

    for(int i=0;i<nn;++i){
      double mi=max(1e-12, vars["m"][i]);
      vars["vvx"][i]=vars["mom_vx"][i]/mi;
      vars["vvy"][i]=vars["mom_vy"][i]/mi;
    }
  }

  void step2_external_forces(){
    // optional: add gravity or drag here (currently none)
  }

  void step3_update_F_compute_stress_and_internal_forces(){
    // gradients from nodal velocity
    mpm->g2pd(vars, {"vvx"}, {"dvx_dx"}, {"dvx_dy"}, false);
    mpm->g2pd(vars, {"vvy"}, {"dvy_dx"}, {"dvy_dy"}, false);

    // update F and compute Piola stress
    for(Part::idx_t p=0;p<mpm->num_particles;++p){
      double L00=clampd(mpm->dp("dvx_dx",p), -P.Lmax, P.Lmax);
      double L01=clampd(mpm->dp("dvx_dy",p), -P.Lmax, P.Lmax);
      double L10=clampd(mpm->dp("dvy_dx",p), -P.Lmax, P.Lmax);
      double L11=clampd(mpm->dp("dvy_dy",p), -P.Lmax, P.Lmax);

      double F00=mpm->dp("F00",p), F01=mpm->dp("F01",p);
      double F10=mpm->dp("F10",p), F11=mpm->dp("F11",p);

      double nF00=(1.0 + P.dt*L00)*F00 + (P.dt*L01)*F10;
      double nF01=(1.0 + P.dt*L00)*F01 + (P.dt*L01)*F11;
      double nF10=(P.dt*L10)*F00 + (1.0 + P.dt*L11)*F10;
      double nF11=(P.dt*L10)*F01 + (1.0 + P.dt*L11)*F11;

      double Jraw=det2(nF00,nF01,nF10,nF11);
      if(!isfinite(Jraw) || Jraw<=1e-12){
        nF00=F00; nF01=F01; nF10=F10; nF11=F11;
        Jraw=det2(nF00,nF01,nF10,nF11);
      }

      double Jtgt=clampd(Jraw,P.Jmin,P.Jmax);

      double inv_sqrtJ=1.0/sqrt(max(1e-12,Jraw));
      nF00*=inv_sqrtJ; nF01*=inv_sqrtJ;
      nF10*=inv_sqrtJ; nF11*=inv_sqrtJ;

      double s=sqrt(max(1e-12,Jtgt));
      nF00*=s; nF01*=s;
      nF10*=s; nF11*=s;

      double J=det2(nF00,nF01,nF10,nF11);
      if(!isfinite(J) || J<=1e-12){
        nF00=1.0; nF01=0.0; nF10=0.0; nF11=1.0;
        J=1.0;
      }

      mpm->dp("F00",p)=nF00; mpm->dp("F01",p)=nF01;
      mpm->dp("F10",p)=nF10; mpm->dp("F11",p)=nF11;
      mpm->dp("J",p)=J;

      double iT00,iT01,iT10,iT11;
      invT2(nF00,nF01,nF10,nF11,iT00,iT01,iT10,iT11);

      double phi_p=clampd(mpm->dp("phi",p),0.0,1.0);
      double soften=pow(max(0.0,1.0-phi_p),P.soften_p);

      double mu_eff=max(P.mu_min, P.mu_s*soften);
      double K0=(P.bulk_ratio*P.mu_s);
      double K_eff=max(P.K_min, K0*soften);

      double a=mu_eff;
      double b=K_eff*(J-1.0)*J;

      mpm->dp("P00",p)=a*(nF00-iT00)+b*iT00;
      mpm->dp("P01",p)=a*(nF01-iT01)+b*iT01;
      mpm->dp("P10",p)=a*(nF10-iT10)+b*iT10;
      mpm->dp("P11",p)=a*(nF11-iT11)+b*iT11;
    }

    // internal force projection
    for(Part::idx_t p=0;p<mpm->num_particles;++p){
      double Vp=mpm->dp("V0",p)*mpm->dp("J",p);
      mpm->dp("Vp_scale",p)=-Vp;
    }

    mpm->p2gd(vars, {"P00"}, {"P01"}, "Vp_scale", {"F_int_vx"}, false);
    mpm->p2gd(vars, {"P10"}, {"P11"}, "Vp_scale", {"F_int_vy"}, false);
  }

  void step4_momentum_balance_compute_nodal_va(){
    for(int i=0;i<nn;++i){
      vars["Ftot_vx"][i]=vars["F_ext_vx"][i]+vars["F_int_vx"][i];
      vars["Ftot_vy"][i]=vars["F_ext_vy"][i]+vars["F_int_vy"][i];

      vars["mom_vx"][i]+=P.dt*vars["Ftot_vx"][i];
      vars["mom_vy"][i]+=P.dt*vars["Ftot_vy"][i];

      double mi=vars["m"][i];
      if(mi>1e-12){
        vars["vvx"][i]=vars["mom_vx"][i]/mi;
        vars["vvy"][i]=vars["mom_vy"][i]/mi;
        vars["avx"][i]=vars["Ftot_vx"][i]/mi;
        vars["avy"][i]=vars["Ftot_vy"][i]/mi;
      }else{
        vars["vvx"][i]=vars["vvy"][i]=0.0;
        vars["avx"][i]=vars["avy"][i]=0.0;
      }

      if(P.node_damp){
        vars["vvx"][i]*=P.damp;
        vars["vvy"][i]*=P.damp;
      }
    }
  }

  void step5_boundary_conditions(){
    // optional: implement wall BC like code2 (edges)
  }

  void step6_g2p_update_particles(){
    for(Part::idx_t p=0;p<mpm->num_particles;++p){
      mpm->dp("vp_x_old",p)=mpm->dp("vp_x",p);
      mpm->dp("vp_y_old",p)=mpm->dp("vp_y",p);
    }

    mpm->g2p(vars, {"vvx","vvy"}, {"vp_x","vp_y"}, false);

    double alpha=clampd(P.pic_alpha,0.0,1.0);
    for(Part::idx_t p=0;p<mpm->num_particles;++p){
      double vpx_pic=mpm->dp("vp_x",p), vpy_pic=mpm->dp("vp_y",p);
      double vpx_old=mpm->dp("vp_x_old",p), vpy_old=mpm->dp("vp_y_old",p);
      mpm->dp("vp_x",p)=(1.0-alpha)*vpx_old + alpha*vpx_pic;
      mpm->dp("vp_y",p)=(1.0-alpha)*vpy_old + alpha*vpy_pic;
    }

    for(Part::idx_t p=0;p<mpm->num_particles;++p){
      mpm->x[p]+=P.dt*mpm->dp("vp_x",p);
      mpm->y[p]+=P.dt*mpm->dp("vp_y",p);
    }
  }

  void mpm_mechanics_step(){
    step0_connectivity();
    clear_nodal();
    step1_p2g_mass_momentum();
    step2_external_forces();
    step3_update_F_compute_stress_and_internal_forces();
    step4_momentum_balance_compute_nodal_va();
    step5_boundary_conditions();
    step6_g2p_update_particles();
  }

  // --------------------------------------------------------------------------
  // Phase-field + heat (replaces enthalpy)
  // --------------------------------------------------------------------------
  void phasefield_thermal_step_allen_cahn(){
    const int NX=P.NX, NY=P.NY;
    const int N=(int)grid.size();

    double Urel_tmp=0.0;
    double h_eff=compute_h_eff_from_relative_velocity(Urel_tmp);

    // hot inflow
    if(P.advect && P.inflow_hot){
      for(int j=0;j<NY;++j){
        int id0=I2(0,j,NX);
        if(grid[id0].mask<=0.0){
          grid[id0].T=P.T_inf;
          grid[id0].phi=0.0;
        }
      }
    }

    vector<double> Tn(N,0.0), phin(N,0.0);
    for(int id=0;id<N;++id){ Tn[id]=grid[id].T; phin[id]=grid[id].phi; }

    // convection source
    vector<double> Qconv_vol(N,0.0);
    if(P.use_conv){
      const double Ax=P.dy*P.thickness;
      const double Ay=P.dx*P.thickness;
      const double V =P.dx*P.dy*P.thickness;

      auto isGas=[&](const Cell& cc){ return cc.mask<=0.0; };
      auto isSolid=[&](const Cell& cc){ return cc.mask>0.0; };

      auto add_pair=[&](int idG,int idS,double Aface){
        double Tg=Tn[idG], Ts=Tn[idS];
        double q=h_eff*(Tg-Ts)*Aface;
        Qconv_vol[idG]+=(-q)/V;
        Qconv_vol[idS]+=(+q)/V;
      };

      for(int j=0;j<NY;++j)
        for(int i=0;i<NX-1;++i){
          int idL=I2(i,j,NX), idR=I2(i+1,j,NX);
          if(isGas(grid[idL]) && isSolid(grid[idR])) add_pair(idL,idR,Ax);
          if(isSolid(grid[idL]) && isGas(grid[idR])) add_pair(idR,idL,Ax);
        }

      for(int j=0;j<NY-1;++j)
        for(int i=0;i<NX;++i){
          int idD=I2(i,j,NX), idU=I2(i,j+1,NX);
          if(isGas(grid[idD]) && isSolid(grid[idU])) add_pair(idD,idU,Ay);
          if(isSolid(grid[idD]) && isGas(grid[idU])) add_pair(idU,idD,Ay);
        }

      for(int id=0;id<N;++id)
        Qconv_vol[id]=clampd(Qconv_vol[id], -P.qconv_max_Wm3, P.qconv_max_Wm3);
    }

    auto k_face=[&](double kL,double kR){ return 2.0*kL*kR/max(1e-12,(kL+kR)); };
    const double inv_dx2=1.0/(P.dx*P.dx);
    const double inv_dy2=1.0/(P.dy*P.dy);

    auto idc=[&](int ii,int jj){
      ii=clampi(ii,0,NX-1); jj=clampi(jj,0,NY-1);
      return I2(ii,jj,NX);
    };

    auto Tnbr=[&](int nid)->double{
      if(P.dirichlet_gas && grid[nid].mask<=0.0) return P.T_inf;
      return Tn[nid];
    };

    // (1) Allen–Cahn
    vector<double> dphi_dt(N,0.0), phinew(N,0.0);

    for(int j=0;j<NY;++j){
      for(int i=0;i<NX;++i){
        int id=I2(i,j,NX);

        if(grid[id].mask<=0.0){
          dphi_dt[id]=0.0;
          phinew[id]=0.0;
          continue;
        }

        double phiC=phin[id];
        double lap =
          (phin[idc(i+1,j)] - 2.0*phiC + phin[idc(i-1,j)])*inv_dx2 +
          (phin[idc(i,j+1)] - 2.0*phiC + phin[idc(i,j-1)])*inv_dy2;

        double dF_dphi = Wprime(phiC) - (P.eps_phi*P.eps_phi)*lap
                         + P.lambda_phi*(Tn[id]-P.Tm)*gprime(phiC);

        dphi_dt[id] = -P.M_phi * dF_dphi;
        phinew[id]  = clampd(phiC + P.dt*dphi_dt[id], P.phi_min, P.phi_max);
      }
    }

    // (2) Heat
    vector<double> Tnew(N,0.0);

    for(int j=0;j<NY;++j){
      for(int i=0;i<NX;++i){
        int id=I2(i,j,NX);

        if(P.dirichlet_gas && grid[id].mask<=0.0){
          Tnew[id]=P.T_inf;
          continue;
        }

        const Cell& c=grid[id];

        double diff=0.0;
        if(i>0){
          int il=I2(i-1,j,NX);
          diff += k_face(c.k,grid[il].k)*(Tnbr(il)-Tn[id])*inv_dx2;
        }
        if(i<NX-1){
          int ir=I2(i+1,j,NX);
          diff += k_face(c.k,grid[ir].k)*(Tnbr(ir)-Tn[id])*inv_dx2;
        }
        if(j>0){
          int jd=I2(i,j-1,NX);
          diff += k_face(c.k,grid[jd].k)*(Tnbr(jd)-Tn[id])*inv_dy2;
        }
        if(j<NY-1){
          int ju=I2(i,j+1,NX);
          diff += k_face(c.k,grid[ju].k)*(Tnbr(ju)-Tn[id])*inv_dy2;
        }

        // advection in gas (simple upwind on T)
        double adv_div=0.0;
        if(P.advect){
          Vec3 u=cell_velocity(c);
          auto T_at=[&](int ii,int jj)->double{ return Tn[idc(ii,jj)]; };

          double Fx_m=(u.x>=0)? u.x*T_at(i-1,j) : u.x*T_at(i,j);
          double Fx_p=(u.x>=0)? u.x*T_at(i,j)   : u.x*T_at(i+1,j);
          double Fy_m=(u.y>=0)? u.y*T_at(i,j-1) : u.y*T_at(i,j);
          double Fy_p=(u.y>=0)? u.y*T_at(i,j)   : u.y*T_at(i,j+1);

          adv_div=(Fx_p-Fx_m)/P.dx + (Fy_p-Fy_m)/P.dy;
        }

        double latent=0.0;
        if(c.mask>0.0){
          latent = P.rho_cr * P.Lm * dphi_dt[id]; // [W/m^3]
        }

        double Qconv = P.use_conv ? Qconv_vol[id] : 0.0;

        double rhoCp=max(1e-12,c.rhoCp);
        double dTdt=(diff - rhoCp*adv_div + latent + Qconv)/rhoCp;

        Tnew[id]=clampd(Tn[id] + P.dt*dTdt, 0.0, P.T_inf);
      }
    }

    // write back
    for(int id=0;id<N;++id){
      grid[id].phi = (grid[id].mask>0.0) ? phinew[id] : 0.0;
      grid[id].T   = Tnew[id];
      update_cell_properties(grid[id]);
    }

    // sample to particles
    for(Part::idx_t p=0;p<mpm->num_particles;++p){
      int i=clampi((int)floor(mpm->x[p]/P.dx),0,NX-1);
      int j=clampi((int)floor(mpm->y[p]/P.dy),0,NY-1);
      const Cell& C=grid[I2(i,j,NX)];
      mpm->dp("T",p)=C.T;
      mpm->dp("phi",p)=C.phi;
    }
  }

  // --------------------------------------------------------------------------
  // VTK output (same as code1)
  // --------------------------------------------------------------------------
  void write_vti_grid_2d(int s){
    string fn=P.out_prefix+"_grid_"+pad5(s)+".vti";
    ofstream f(fn);
    if(!f){ cerr<<"Cannot open "<<fn<<"\n"; return; }

    f<<"<?xml version=\"1.0\"?>\n";
    f<<"<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    f<<"  <ImageData Origin=\"0 0 0\" Spacing=\""<<P.dx<<" "<<P.dy<<" "<<P.thickness
     <<"\" WholeExtent=\"0 "<<(P.NX-1)<<" 0 "<<(P.NY-1)<<" 0 0\">\n";
    f<<"    <Piece Extent=\"0 "<<(P.NX-1)<<" 0 "<<(P.NY-1)<<" 0 0\">\n";
    f<<"      <CellData Scalars=\"T\">\n";

    auto write_scalar=[&](const string& name, auto getter){
      f<<"        <DataArray type=\"Float32\" Name=\""<<name
       <<"\" format=\"ascii\" NumberOfComponents=\"1\">\n";
      for(int j=0;j<P.NY;++j)
        for(int i=0;i<P.NX;++i)
          f<<(float)getter(grid[I2(i,j,P.NX)])<<"\n";
      f<<"        </DataArray>\n";
    };

    write_scalar("T",   [](const Cell& c){return c.T;});
    write_scalar("phi", [](const Cell& c){return c.phi;});
    write_scalar("mask",[](const Cell& c){return c.mask;});

    f<<"      </CellData>\n";
    f<<"      <PointData></PointData>\n";
    f<<"    </Piece>\n";
    f<<"  </ImageData>\n";
    f<<"</VTKFile>\n";
    cerr<<"Wrote "<<fn<<"\n";
  }

  void write_vtp_particles(int s){
    string fn=P.out_prefix+"_particles_"+pad5(s)+".vtp";
    ofstream f(fn);
    if(!f){ cerr<<"Cannot open "<<fn<<"\n"; return; }

    const size_t Np=(size_t)mpm->num_particles;

    f<<"<?xml version=\"1.0\"?>\n";
    f<<"<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    f<<"  <PolyData>\n";
    f<<"    <Piece NumberOfPoints=\""<<Np<<"\" NumberOfVerts=\""<<Np
     <<"\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n";

    f<<"      <Points>\n";
    f<<"        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for(size_t p=0;p<Np;++p) f<<(float)mpm->x[p]<<" "<<(float)mpm->y[p]<<" "<<0.0f<<"\n";
    f<<"        </DataArray>\n";
    f<<"      </Points>\n";

    f<<"      <Verts>\n";
    f<<"        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for(size_t i=0;i<Np;++i) f<<(int)i<<"\n";
    f<<"        </DataArray>\n";
    f<<"        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for(size_t i=0;i<Np;++i) f<<(int)(i+1)<<"\n";
    f<<"        </DataArray>\n";
    f<<"      </Verts>\n";

    f<<"      <PointData Scalars=\"T\">\n";
    auto write_pscalar=[&](const string& name, auto getter){
      f<<"        <DataArray type=\"Float32\" Name=\""<<name
       <<"\" format=\"ascii\" NumberOfComponents=\"1\">\n";
      for(size_t p=0;p<Np;++p) f<<(float)getter((Part::idx_t)p)<<"\n";
      f<<"        </DataArray>\n";
    };

    write_pscalar("T",[&](Part::idx_t p){return mpm->dp("T",p);});
    write_pscalar("phi",[&](Part::idx_t p){return mpm->dp("phi",p);});
    write_pscalar("J",[&](Part::idx_t p){return mpm->dp("J",p);});

    f<<"      </PointData>\n";
    f<<"      <CellData></CellData>\n";
    f<<"    </Piece>\n";
    f<<"  </PolyData>\n";
    f<<"</VTKFile>\n";
    cerr<<"Wrote "<<fn<<"\n";
  }

  // --------------------------------------------------------------------------
  // run
  // --------------------------------------------------------------------------
  void run(){
    init_thermal_fields_once();

    for(int s=0;s<P.steps;++s){
      mpm_mechanics_step();

      if(s % P.rebuild_every == 0)
        rebuild_mask_phi_T_from_particles();

      phasefield_thermal_step_allen_cahn();

      if(s % P.output_every == 0){
        double phim=0.0, Tavg=0.0, Javg=0.0;
        for(Part::idx_t p=0;p<mpm->num_particles;++p){
          phim += mpm->dp("phi",p);
          Tavg += mpm->dp("T",p);
          Javg += mpm->dp("J",p);
        }
        if(mpm->num_particles>0){
          phim /= (double)mpm->num_particles;
          Tavg /= (double)mpm->num_particles;
          Javg /= (double)mpm->num_particles;
        }

        double Urel_dbg=0.0;
        double h_eff_dbg=compute_h_eff_from_relative_velocity(Urel_dbg);

        cout<<"step "<<s
            <<" phi_avg="<<phim
            <<" T_avg="<<Tavg
            <<" J_avg="<<Javg
            <<" Urel="<<Urel_dbg
            <<" h_eff="<<h_eff_dbg
            <<"\n";

        write_vti_grid_2d(s);
        write_vtp_particles(s);
      }
    }
  }
};

// --------------------------------------------------------------------------
// Map DATA -> Params (like code2 style). Adjust if your DATA fields differ.
// --------------------------------------------------------------------------
static Params params_from_DATA(const DATA& data){
  Params P;

  P.NX = data.NX;
  P.NY = data.NY;
  P.dx = data.dx;
  P.dy = data.dy;
  P.thickness = data.thickness;

  P.dt = data.dt;
  P.steps = data.steps;
  P.output_every = data.output_every;
  P.rebuild_every = data.rebuild_every;

  P.T_inf = data.T_inf;
  P.u_gas = Vec3{data.u_gas.x, data.u_gas.y, data.u_gas.z};
  P.advect = data.advect;
  P.inflow_hot = data.inflow_hot;
  P.dirichlet_gas = data.dirichlet_gas;

  P.use_initial_shear = data.use_initial_shear;
  P.shear_v0 = data.shear_v0;

  P.M_phi = data.M_phi;
  P.eps_phi = data.eps_phi;
  P.lambda_phi = data.lambda_phi;

  P.out_prefix = data.out_prefix;

  return P;
}


int main(int argc, char** argv){
  string json_path = (argc>1) ? argv[1] : string("DATA.json");

  DATA data(json_path.c_str());   // like your code2
  Params P = params_from_DATA(data);

  Sim sim(P);
sim.seed_particles_disk_from_DATA(data);
sim.init_thermal_fields_once();
sim.run();

  return 0;
}

