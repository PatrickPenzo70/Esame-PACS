// ============================================================================
// 2D MPM + Enthalpy Melting Demo (Cr2O3 disk in hot Ar+H2 gas)
// CASE A) DEFORMATION ENABLED (initial shear) + THERMO FIX (correct P2G(T))
// ============================================================================

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "quadgrid_cpp.h"
#include "particles.h"

using namespace std;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

struct Vec3 {
  double x=0,y=0,z=0;
  Vec3(){}
  Vec3(double X,double Y,double Z):x(X),y(Y),z(Z){}
};

static inline int clampi(int a,int lo,int hi){ return max(lo,min(a,hi)); }
static inline double clampd(double x,double lo,double hi){ return max(lo,min(x,hi)); }
static inline double det2(double a00,double a01,double a10,double a11){ return a00*a11 - a01*a10; }

static inline void invT2(double F00,double F01,double F10,double F11,
                         double &iT00,double &iT01,double &iT10,double &iT11){
  double J=det2(F00,F01,F10,F11);
  double invJ=1.0/max(1e-12,J);
  iT00= F11*invJ;  iT01=-F01*invJ;
  iT10=-F10*invJ;  iT11= F00*invJ;
}

static inline int I2(int i,int j,int NX){ return j*NX+i; }

struct Params {
  int NX=192, NY=192;
  double dx=1e-5, dy=1e-5;
  double thickness=1e-5;
  double dt=2e-10;
  int steps=200000;
  int output_every=200;

  // Gas
  double rho_g=0.2;
  double cp_g=1500.0;
  double k_g=1.0;
  double T_inf=15000.0;
  Vec3 u_gas{300.0,0.0,0.0};
  bool advect=true;
  bool inflow_hot=false;     // off by default
  bool dirichlet_gas=false;  // off by default

  // Solid (Cr2O3)
  double rho_cr=5200.0;
  double cp_s=800.0;
  double cp_l=1000.0;
  double k_s=5.0;
  double k_l=3.0;
  double Tm=2435.0;
  double Lm=4e5;
  
  // Coupling melt->mechanics
  double soften_p = 3.0;     // esponente p
  double mu_min = 1e4;       // Pa (quasi fluido, non zero per stabilità)
  double K_min  = 1e5;       // Pa

  Vec3 center{0.00096,0.00096,0.0};
  double R=100e-6;
  int nParticles=15000;
  
  double pic_alpha = 0.05; // 0 = FLIP puro (instabile), 1 = PIC puro (molto diffuso)

  // gas<->solid convection (reduced)
  bool use_conv=true;
  double h_conv=5e3;
  double qconv_max_Wm3=5e11;

  // MPM mechanics
  double mu_s=5e8;
  double bulk_ratio=200.0;
  double Jmin=0.60;
  double Jmax=1.40;

  bool node_damp=true;
  double damp=0.99;
  double Lmax=2e6;

  // CASE A initial shear
  bool use_initial_shear=true;
  double shear_v0=2.0;
};

struct Cell{
  double T=300.0;
  double rhoCp=0.0;
  double k=0.0;
  double phi=0.0;
  double mask=0.0;
};

struct Sim{
  using QGrid = quadgrid_t<std::vector<double>>;
  using Part  = particles_t;

  Params P;
  vector<Cell> grid;

  QGrid qg;
  unique_ptr<Part> mpm;

  map<string, vector<double>> gv;
  int nn=0;

  Sim(const Params& p):P(p){
    grid.resize((size_t)P.NX*P.NY);
    qg.set_sizes((QGrid::idx_t)P.NY,(QGrid::idx_t)P.NX,P.dy,P.dx);
    nn=(int)qg.num_local_nodes();

    gv["m"]=vector<double>(nn,0.0);
    gv["vx"]=vector<double>(nn,0.0);
    gv["vy"]=vector<double>(nn,0.0);
    gv["fx"]=vector<double>(nn,0.0);
    gv["fy"]=vector<double>(nn,0.0);
  }

  static string pad5(int s){
    ostringstream oss; oss<<setw(5)<<setfill('0')<<s; return oss.str();
  }

  Vec3 cell_velocity(const Cell& c) const {
    if(c.mask>0.0) return Vec3(0,0,0);
    return P.u_gas;
  }

  void update_cell_properties(Cell& c){
    if(c.mask>0.0){
      double f=c.mask;
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

  void seed_particles_quadgrid(){
    mt19937 rng(42);
    uniform_real_distribution<double> U(0.0,1.0);

    vector<double> xv(P.nParticles), yv(P.nParticles);

    double Adisk=M_PI*P.R*P.R;
    double Apt=Adisk/(double)P.nParticles;
    double Vpt=Apt*P.thickness;
    double mpt=P.rho_cr*Vpt;

    for(int n=0;n<P.nParticles;++n){
      double r=P.R*sqrt(U(rng));
      double a=2.0*M_PI*U(rng);
      xv[n]=P.center.x + r*cos(a);
      yv[n]=P.center.y + r*sin(a);
    }

    vector<string> iprops={};
    vector<string> dprops={
      "mp","vp_x","vp_y","V0",
      "F00","F01","F10","F11","J",
      "T","phi",
      "dvx_dx","dvx_dy","dvy_dx","dvy_dy",
      "P00","P01","P10","P11"
    };

    mpm=make_unique<Part>((Part::idx_t)P.nParticles, iprops, dprops, qg, xv, yv);
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
        mpm->dp("vp_x",p)= P.shear_v0*((yy-P.center.y)/max(1e-12,P.R));
        mpm->dp("vp_y",p)= 0.0;
      }else{
        mpm->dp("vp_x",p)=0.0;
        mpm->dp("vp_y",p)=0.0;
      }
    }
  }

  // ------------------ IMPORTANT FIX: always P2G(T,phi) for solid cells ------------------

  void init_thermal_fields_once(){
    for(auto& c: grid){
      c.mask=0.0; c.phi=0.0;
      c.T=P.T_inf;                 // gas hot initial condition
      c.rhoCp=P.rho_g*P.cp_g;
      c.k=P.k_g;
    }
    rebuild_mask_phi_T_from_particles(); // put solid T=300 initially
  }

  void rebuild_mask_phi_T_from_particles(){
    const double Vcell=P.dx*P.dy*P.thickness;

    vector<double> mcell((size_t)P.NX*P.NY,0.0);
    vector<double> mp_sum((size_t)P.NX*P.NY,0.0);
    vector<double> mpT_sum((size_t)P.NX*P.NY,0.0);
    vector<double> mpPhi_sum((size_t)P.NX*P.NY,0.0);

    // reset mask/phi only (do not reset T for gas)
    for(auto& c: grid){ c.mask=0.0; c.phi=0.0; }

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
          // solid cell: take temperature + phi from particles (mass-weighted)
          c.T   = mpT_sum[id]/mp_sum[id];
          c.phi = clampd(mpPhi_sum[id]/mp_sum[id],0.0,1.0);
        }else{
          // gas cell: keep its current T (do not overwrite here)
          c.phi = 0.0;
          if(P.dirichlet_gas) c.T = P.T_inf;
        }

        update_cell_properties(c);
      }
    }
  }

  // ------------------ MPM hyperelastic step (same as before) ------------------

  void mpm_step_hyperelastic(){
    fill(gv["m"].begin(), gv["m"].end(), 0.0);
    fill(gv["vx"].begin(), gv["vx"].end(), 0.0);
    fill(gv["vy"].begin(), gv["vy"].end(), 0.0);
    fill(gv["fx"].begin(), gv["fx"].end(), 0.0);
    fill(gv["fy"].begin(), gv["fy"].end(), 0.0);

    mpm->p2g(gv, {"mp"}, {"m"}, false);

    if(mpm->dprops.find("px")==mpm->dprops.end()){
      mpm->dprops["px"].resize(mpm->num_particles);
      mpm->dprops["py"].resize(mpm->num_particles);
    }
    for(Part::idx_t p=0;p<mpm->num_particles;++p){
      mpm->dp("px",p)=mpm->dp("mp",p)*mpm->dp("vp_x",p);
      mpm->dp("py",p)=mpm->dp("mp",p)*mpm->dp("vp_y",p);
    }
    mpm->p2g(gv, {"px","py"}, {"vx","vy"}, false);

    for(int i=0;i<nn;++i){
      double mi=max(1e-12, gv["m"][i]);
      gv["vx"][i]/=mi;
      gv["vy"][i]/=mi;
    }

    mpm->g2pd(gv, {"vx"}, {"dvx_dx"}, {"dvx_dy"}, false);
    mpm->g2pd(gv, {"vy"}, {"dvy_dx"}, {"dvy_dy"}, false);

    auto mu_of_phi = [&](double phi){
    double a = pow(max(0.0, 1.0 - phi), P.soften_p);
    return max(P.mu_min, P.mu_s * a);
    };
    auto K_of_phi = [&](double phi, double mu_eff){
    double Ks = P.bulk_ratio * P.mu_s;
    double a  = pow(max(0.0, 1.0 - phi), P.soften_p);
    double K_eff = Ks * a;
    return max(P.K_min, K_eff);
    };


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

      double Jtgt = clampd(Jraw, P.Jmin, P.Jmax);

      // se Jraw è già fuori range, NON scalare soltanto: resetta volumetria in modo più robusto
      // 1) rendi isocoro: F <- F / sqrt(Jraw)  (2D)
      double inv_sqrtJ = 1.0 / sqrt(max(1e-12, Jraw));
      nF00 *= inv_sqrtJ; nF01 *= inv_sqrtJ;
      nF10 *= inv_sqrtJ; nF11 *= inv_sqrtJ;

      // ora J=1 (circa). re-imponi Jtgt moltiplicando
      double s = sqrt(max(1e-12, Jtgt));
      nF00 *= s; nF01 *= s;
      nF10 *= s; nF11 *= s;

      // safety finale: se ancora non va, torna a identità
      double Jchk = det2(nF00,nF01,nF10,nF11);
      if(!isfinite(Jchk) || Jchk<=1e-12){
        nF00=1.0; nF01=0.0; nF10=0.0; nF11=1.0;
      }


      double J=det2(nF00,nF01,nF10,nF11);
      mpm->dp("F00",p)=nF00; mpm->dp("F01",p)=nF01;
      mpm->dp("F10",p)=nF10; mpm->dp("F11",p)=nF11;
      mpm->dp("J",p)=J;

      double iT00,iT01,iT10,iT11;
      invT2(nF00,nF01,nF10,nF11,iT00,iT01,iT10,iT11);
      
      // === melt -> mechanics coupling (QUI) ===
      double phi_p = clampd(mpm->dp("phi", p), 0.0, 1.0);
      double soften = pow(max(0.0, 1.0 - phi_p), P.soften_p);

      double mu_eff = max(P.mu_min, P.mu_s * soften);
      double K_eff  = max(P.K_min,  (P.bulk_ratio * P.mu_s) * soften);

      // P = mu(F - F^{-T}) + K*(J-1)*J*F^{-T}
      double a = mu_eff;
      double b = K_eff * (J - 1.0) * J;

      mpm->dp("P00",p)=a*(nF00-iT00)+b*iT00;
      mpm->dp("P01",p)=a*(nF01-iT01)+b*iT01;
      mpm->dp("P10",p)=a*(nF10-iT10)+b*iT10;
      mpm->dp("P11",p)=a*(nF11-iT11)+b*iT11;
    }

    if(mpm->dprops.find("Vp_scale")==mpm->dprops.end())
      mpm->dprops["Vp_scale"].resize(mpm->num_particles);

    for(Part::idx_t p=0;p<mpm->num_particles;++p){
      double Vp=mpm->dp("V0",p)*mpm->dp("J",p);
      mpm->dp("Vp_scale",p)=-Vp;
    }

    mpm->p2gd(gv, {"P00"}, {"P01"}, "Vp_scale", {"fx"}, false);
    mpm->p2gd(gv, {"P10"}, {"P11"}, "Vp_scale", {"fy"}, false);

    for(int i=0;i<nn;++i){
      double mi=gv["m"][i];
      if(mi<=0.0) continue;
      gv["vx"][i]+=P.dt*(gv["fx"][i]/mi);
      gv["vy"][i]+=P.dt*(gv["fy"][i]/mi);
      if(P.node_damp){ gv["vx"][i]*=P.damp; gv["vy"][i]*=P.damp; }
    }

    if(mpm->dprops.find("vp_x_old")==mpm->dprops.end()){
     mpm->dprops["vp_x_old"].resize(mpm->num_particles);
     mpm->dprops["vp_y_old"].resize(mpm->num_particles);
    }
    for(Part::idx_t p=0;p<mpm->num_particles;++p){
     mpm->dp("vp_x_old",p)=mpm->dp("vp_x",p);
     mpm->dp("vp_y_old",p)=mpm->dp("vp_y",p);
    }

    mpm->g2p(gv, {"vx","vy"}, {"vp_x","vp_y"}, false);
    
    double a = clampd(P.pic_alpha, 0.0, 1.0);
     for(Part::idx_t p=0;p<mpm->num_particles;++p){
     double vpx_pic = mpm->dp("vp_x",p);
     double vpy_pic = mpm->dp("vp_y",p);
     double vpx_old = mpm->dp("vp_x_old",p);
     double vpy_old = mpm->dp("vp_y_old",p);
    // (1-a)*old + a*PIC  => smorza instabilità
       mpm->dp("vp_x",p) = (1.0-a)*vpx_old + a*vpx_pic;
       mpm->dp("vp_y",p) = (1.0-a)*vpy_old + a*vpy_pic;
     }
 
    for(Part::idx_t p=0;p<mpm->num_particles;++p){
      mpm->x[p]+=P.dt*mpm->dp("vp_x",p);
      mpm->y[p]+=P.dt*mpm->dp("vp_y",p);
    }

    mpm->init_particle_mesh();
  }

  // ------------------ thermo enthalpy + melt (same as your last) ------------------

  void diffuse_and_melt_step_enthalpy_2d(){
    const int NX=P.NX, NY=P.NY;
    const int N=(int)grid.size();

    if(P.advect && P.inflow_hot){
      for(int j=0;j<NY;++j){
        int id0=I2(0,j,NX);
        if(grid[id0].mask<=0.0){
          grid[id0].T=P.T_inf;
          grid[id0].phi=0.0;
        }
      }
    }

    vector<double> H(N,0.0), Hnew(N,0.0);
    for(int id=0;id<N;++id){
      const Cell& c=grid[id];
      double rhoL_eff=(c.mask>0.0)?(P.rho_cr*c.mask*P.Lm):0.0;
      H[id]=c.rhoCp*c.T + rhoL_eff*c.phi;
    }

    vector<double> Qconv_vol(N,0.0);
    if(P.use_conv){
      const double Ax=P.dy*P.thickness;
      const double Ay=P.dx*P.thickness;
      const double V =P.dx*P.dy*P.thickness;

      auto isGas=[&](const Cell& cc){ return cc.mask<=0.0; };
      auto isSolid=[&](const Cell& cc){ return cc.mask>0.0; };

      auto add_pair=[&](int idG,int idS,double Aface){
        double Tg=grid[idG].T, Ts=grid[idS].T;
        double q=P.h_conv*(Tg-Ts)*Aface;
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

    auto Tnbr=[&](int nid)->double{
      if(P.dirichlet_gas && grid[nid].mask<=0.0) return P.T_inf;
      return grid[nid].T;
    };

    for(int j=0;j<NY;++j){
      for(int i=0;i<NX;++i){
        int id=I2(i,j,NX);
        Cell& c=grid[id];

        if(P.dirichlet_gas && c.mask<=0.0){
          Hnew[id]=c.rhoCp*P.T_inf;
          continue;
        }

        double diff=0.0;
        if(i>0){
          int il=I2(i-1,j,NX);
          diff+=k_face(c.k,grid[il].k)*(Tnbr(il)-c.T)*inv_dx2;
        }
        if(i<NX-1){
          int ir=I2(i+1,j,NX);
          diff+=k_face(c.k,grid[ir].k)*(Tnbr(ir)-c.T)*inv_dx2;
        }
        if(j>0){
          int jd=I2(i,j-1,NX);
          diff+=k_face(c.k,grid[jd].k)*(Tnbr(jd)-c.T)*inv_dy2;
        }
        if(j<NY-1){
          int ju=I2(i,j+1,NX);
          diff+=k_face(c.k,grid[ju].k)*(Tnbr(ju)-c.T)*inv_dy2;
        }

        double adv_div=0.0;
        if(P.advect){
          Vec3 u=cell_velocity(c);
          auto H_at=[&](int ii,int jj)->double{
            ii=clampi(ii,0,NX-1); jj=clampi(jj,0,NY-1);
            return H[I2(ii,jj,NX)];
          };

          double Fx_m=0,Fx_p=0,Fy_m=0,Fy_p=0;
          Fx_m=(u.x>=0)? u.x*H_at(i-1,j) : u.x*H_at(i,j);
          Fx_p=(u.x>=0)? u.x*H_at(i,j)   : u.x*H_at(i+1,j);
          Fy_m=(u.y>=0)? u.y*H_at(i,j-1) : u.y*H_at(i,j);
          Fy_p=(u.y>=0)? u.y*H_at(i,j)   : u.y*H_at(i,j+1);

          adv_div=(Fx_p-Fx_m)/P.dx + (Fy_p-Fy_m)/P.dy;
        }

        double Qconv=P.use_conv?Qconv_vol[id]:0.0;
        Hnew[id]=H[id] + P.dt*(diff - adv_div + Qconv);
        if(!isfinite(Hnew[id]) || Hnew[id]<0.0) Hnew[id]=H[id];
      }
    }

    for(int id=0;id<N;++id){
      const Cell& c=grid[id];
      if(c.mask<=0.0){
        Hnew[id]=clampd(Hnew[id],0.0,c.rhoCp*P.T_inf);
      }else{
        double rhoL_eff=P.rho_cr*c.mask*P.Lm;
        Hnew[id]=clampd(Hnew[id],0.0,c.rhoCp*P.T_inf + rhoL_eff);
      }
    }

    for(int j=0;j<NY;++j){
      for(int i=0;i<NX;++i){
        int id=I2(i,j,NX);
        Cell& c=grid[id];

        if(P.dirichlet_gas && c.mask<=0.0){
          c.T=P.T_inf; c.phi=0.0; update_cell_properties(c); continue;
        }

        if(c.mask<=0.0){
          c.T=Hnew[id]/max(1e-12,c.rhoCp);
          c.T=min(c.T,P.T_inf);
          c.phi=0.0;
          update_cell_properties(c);
          continue;
        }

        double rhoCp=max(1e-12,c.rhoCp);
        double rhoL_eff=max(1e-20,P.rho_cr*c.mask*P.Lm);

        double Hs=rhoCp*P.Tm;
        double Hl=rhoCp*P.Tm + rhoL_eff;
        double Hv=Hnew[id];

        if(Hv<Hs){ c.phi=0.0; c.T=Hv/rhoCp; }
        else if(Hv<=Hl){ c.T=P.Tm; c.phi=clampd((Hv-Hs)/rhoL_eff,0.0,1.0); }
        else{ c.phi=1.0; c.T=(Hv-rhoL_eff)/rhoCp; }

        c.T=min(c.T,P.T_inf);
        update_cell_properties(c);
      }
    }

    for(Part::idx_t p=0;p<mpm->num_particles;++p){
      int i=clampi((int)floor(mpm->x[p]/P.dx),0,NX-1);
      int j=clampi((int)floor(mpm->y[p]/P.dy),0,NY-1);
      const Cell& C=grid[I2(i,j,NX)];
      mpm->dp("T",p)=C.T;
      mpm->dp("phi",p)=C.phi;
    }
  }

  // ------------------ VTK output ------------------

  void write_vti_grid_2d(int s){
    string fn=string("mpm2d_grid_")+pad5(s)+".vti";
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
    string fn=string("mpm2d_particles_")+pad5(s)+".vtp";
    ofstream f(fn);
    if(!f){ cerr<<"Cannot open "<<fn<<"\n"; return; }

    const size_t N=(size_t)mpm->num_particles;

    f<<"<?xml version=\"1.0\"?>\n";
    f<<"<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    f<<"  <PolyData>\n";
    f<<"    <Piece NumberOfPoints=\""<<N<<"\" NumberOfVerts=\""<<N
     <<"\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n";

    f<<"      <Points>\n";
    f<<"        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for(size_t p=0;p<N;++p) f<<(float)mpm->x[p]<<" "<<(float)mpm->y[p]<<" "<<0.0f<<"\n";
    f<<"        </DataArray>\n";
    f<<"      </Points>\n";

    f<<"      <Verts>\n";
    f<<"        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for(size_t i=0;i<N;++i) f<<(int)i<<"\n";
    f<<"        </DataArray>\n";
    f<<"        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for(size_t i=0;i<N;++i) f<<(int)(i+1)<<"\n";
    f<<"        </DataArray>\n";
    f<<"      </Verts>\n";

    f<<"      <PointData Scalars=\"T\">\n";
    auto write_pscalar=[&](const string& name, auto getter){
      f<<"        <DataArray type=\"Float32\" Name=\""<<name
       <<"\" format=\"ascii\" NumberOfComponents=\"1\">\n";
      for(size_t p=0;p<N;++p) f<<(float)getter((Part::idx_t)p)<<"\n";
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

  void run(){
    seed_particles_quadgrid();
    init_thermal_fields_once();

    for(int s=0;s<P.steps;++s){
      mpm_step_hyperelastic();

      // IMPORTANT FIX: update mask/phi/T from particles EVERY step
      rebuild_mask_phi_T_from_particles();

      diffuse_and_melt_step_enthalpy_2d();

      if(s%P.output_every==0){
        double phim=0,Tavg=0,Javg=0;
        for(Part::idx_t p=0;p<mpm->num_particles;++p){
          phim+=mpm->dp("phi",p);
          Tavg+=mpm->dp("T",p);
          Javg+=mpm->dp("J",p);
        }
        if(mpm->num_particles>0){
          phim/= (double)mpm->num_particles;
          Tavg/= (double)mpm->num_particles;
          Javg/= (double)mpm->num_particles;
        }

        cout<<"step "<<s<<" phi_avg="<<phim<<" T_avg="<<Tavg<<" J_avg="<<Javg<<"\n";
        write_vti_grid_2d(s);
        write_vtp_particles(s);
      }
    }
  }
};

int main(){
  Params P;
  Sim sim(P);
  sim.run();
  return 0;
}

