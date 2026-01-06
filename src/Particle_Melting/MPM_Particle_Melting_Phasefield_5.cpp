#include <bits/stdc++.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// -----------------------------------------------------------------------------
// MPM (meccanico) + TERMICA + FUSIONE con quadgrid_cpp.h / particles.h
// - Griglia regolare 2D da quadgrid_t (gestione indici, dx, dy, VTK export)
// - Particelle: x,y,vx,vy, F (matrice 2x2), T, phi
// - Meccanica: P2G/G2P esplicito bilineare (hat), Neo-Hookean semplice, gravità
// - Termica: diffusione su griglia con k variabile, sorgente latente S = -rho Lm dphi
// - Mask (frazione solida per cella) ricostruita da volume particellare
// - OpenMP opzionale; VTK export di campi termici
// -----------------------------------------------------------------------------

#include <quadgrid_cpp.h>
#include <particles.h>
#include <memory>

struct Mat2{ double a11=1, a12=0, a21=0, a22=1; };
static inline Mat2 I2(){ return {1,0,0,1}; }
static inline double det(const Mat2&A){ return A.a11*A.a22 - A.a12*A.a21; }

struct Params{
    // Griglia
    int NX=200, NY=200;      // nodi per lato (quadgrid usa nodi regolari)
    double dx=1e-5, dy=1e-5; // passo [m]

    // Tempo
    int steps=6000;
    int output_every=200;
    double cfl_mech=0.25;    // CFL meccanica

    // Materiale Cr2O3
    double rho=5200.0;          // kg/m^3
    double E=120e9;             // Pa (Young)
    double nu=0.25;             // Poisson

    // Termica
    double k_s=5.0, k_l=3.0;          // W/mK
    double cp_s=800.0, cp_l=1000.0;   // J/kgK
    double Tm=2435.0;                 // K
    double Lm=4e5;                    // J/kg
    double Kphi=1e-3;                 // 1/s/K
    
    // Phase-field (Allen–Cahn)
    double Mphi     = 5e-6;   // mobilità [1/s] (tuning numerico)
    double eps_phi  = 2.5e-5; // spessore interfaccia ~ 2-5*dx [m]
    double W_phi    = 1.0;    // barriera doppio pozzo (tuning)
    double lam_phi  = 5e-4;   // accoppiamento termico (tuning)

    // Gas (ambiente)
    double rho_g=0.2, cp_g=1500.0, k_g=10.0, T_inf=15000.0;
    bool dirichlet_gas=true;

    // Gravità e PIC/FLIP
    double gx=0.0, gy=-9.81;  // m/s^2
    double picflip=0.10;      // 0..1

    // Particella iniziale
    double cx=0.001, cy=0.001; // centro
    double R = 100e-6;         // raggio [m]
    int nParticles = 8000;
};

using grid_t = quadgrid_t<std::vector<double>>;
using idx_t  = grid_t::idx_t;

struct MechGrid{
    // nodi (stessa numerazione di quadgrid)
    std::vector<double> m;     // massa
    std::vector<double> px, py;// quantità di moto
    std::vector<double> fx, fy;// forze
    std::vector<double> vx, vy;// velocità nodali

    explicit MechGrid(size_t N): m(N,0), px(N,0), py(N,0), fx(N,0), fy(N,0), vx(N,0), vy(N,0){}
    void clear(){
        std::fill(m.begin(), m.end(), 0.0);
        std::fill(px.begin(), px.end(), 0.0);
        std::fill(py.begin(), py.end(), 0.0);
        std::fill(fx.begin(), fx.end(), 0.0);
        std::fill(fy.begin(), fy.end(), 0.0);
        std::fill(vx.begin(), vx.end(), 0.0);
        std::fill(vy.begin(), vy.end(), 0.0);
    }
};

struct MPMThermoSim{
    Params P;
    grid_t grid;                                  // griglia per termo + indici
    MechGrid g;                                   // griglia meccanica
    std::unique_ptr<particles_t> pt;              // particelle (x,y,vx,vy,T,phi, ecc.)
    std::vector<Mat2> F;                          // F per particella
    std::vector<double> vol0;                     // volume iniziale per particella (2D: area cella come spessore)
    std::vector<double> massp;                    // massa particella

    // campi termici su griglia (per VTK e calcolo)
    std::map<std::string, std::vector<double>> vars; // {T, phi, mask, rhoCp, k, S}

    MPMThermoSim(const Params& p): P(p), g(0), pt(nullptr){
        // 1) griglia
        grid.set_sizes(P.NY, P.NX, P.dx, P.dy);
        g = MechGrid(grid.num_global_nodes());

        // 2) particelle: disco
        std::vector<double> x(P.nParticles), y(P.nParticles);
        std::mt19937 rng(42); std::uniform_real_distribution<double> U(0,1);
        for(int i=0;i<P.nParticles;++i){
            double r = P.R*std::sqrt(U(rng)); double th = 2*M_PI*U(rng);
            x[i] = P.cx + r*std::cos(th);
            y[i] = P.cy + r*std::sin(th);
        }
        pt = std::make_unique<particles_t>(
          P.nParticles,
          std::vector<std::string>{"id"},
          std::vector<std::string>{"T","phi","vx","vy"},
          grid,
          x, y
        );
        pt.reset(new particles_t(
    P.nParticles,
    std::vector<std::string>{"id"},
    std::vector<std::string>{"T","phi","vx","vy"},
    grid,
    x, y
));

        for (idx_t i=0; i<(idx_t)pt->x.size(); ++i){
            pt->dprops["T"][i]=300.0;
            pt->dprops["phi"][i]=0.0;
            pt->dprops["vx"][i]=0.0;
            pt->dprops["vy"][i]=0.0;
        }

        F.assign(P.nParticles, I2());
        vol0.assign(P.nParticles, P.dx*P.dy); // spessore=1
        massp.assign(P.nParticles, P.rho * vol0[0]);

        // 3) campi termici
        vars = {
            {"T",     std::vector<double>(grid.num_global_nodes(), P.T_inf)},
            {"phi",   std::vector<double>(grid.num_global_nodes(), 0.0)},
            {"mask",  std::vector<double>(grid.num_global_nodes(), 0.0)},
            {"rhoCp", std::vector<double>(grid.num_global_nodes(), 0.0)},
            {"k",     std::vector<double>(grid.num_global_nodes(), 0.0)},
            {"S",     std::vector<double>(grid.num_global_nodes(), 0.0)},
            {"phidot",std::vector<double>(grid.num_global_nodes(), 0.0)} 
        };

        // inizializza maschera/termica da particelle
        rebuild_mask_from_particles();
        set_thermal_properties();
        // set T fredda dove c'è particella
        for (idx_t id=0; id<grid.num_global_nodes(); ++id)
            if (vars["mask"][id]>0) vars["T"][id]=300.0;
    }

    inline idx_t I(int i,int j) const { return j*P.NX + i; }

    // --------- Meccanica ---------
    void p2g_mech(){
        const double mu = P.E/(2*(1+P.nu));
        const double la = P.E*P.nu/((1+P.nu)*(1-2*P.nu));
        g.clear();

        #pragma omp parallel for
        for(long long p=0; p<(long long)(idx_t)pt->x.size(); ++p){
            // cella base
            double xp = pt->x[p], yp = pt->y[p];
            int i0 = (int)std::floor(xp / P.dx);
            int j0 = (int)std::floor(yp / P.dy);
            double fx = xp / P.dx - i0; if(i0<0){i0=0; fx=0;} if(i0>P.NX-2){i0=P.NX-2; fx=1;}
            double fy = yp / P.dy - j0; if(j0<0){j0=0; fy=0;} if(j0>P.NY-2){j0=P.NY-2; fy=1;}
            double wx[2] = {1.0 - fx, fx};
            double wy[2] = {1.0 - fy, fy};
            double dwx[2] = {-1.0/P.dx, 1.0/P.dx};
            double dwy[2] = {-1.0/P.dy, 1.0/P.dy};

            // stress Neo-Hookean: P = mu(F - F^{-T}) + la lnJ F^{-T}
            Mat2 Fp = F[p];
            double J = det(Fp);
            double invT11 =  Fp.a22/J, invT12 = -Fp.a21/J;
            double invT21 = -Fp.a12/J, invT22 =  Fp.a11/J;
            // 1° Piola
            double P11 = mu*(Fp.a11 - invT11) + la*std::log(J)*invT11;
            double P12 = mu*(Fp.a12 - invT12) + la*std::log(J)*invT12;
            double P21 = mu*(Fp.a21 - invT21) + la*std::log(J)*invT21;
            double P22 = mu*(Fp.a22 - invT22) + la*std::log(J)*invT22;

            double m = massp[p];
            double vx = pt->dprops["vx"][p];
            double vy = pt->dprops["vy"][p];
            double Vp = vol0[p];

            for(int a=0;a<2;++a){
                for(int b=0;b<2;++b){
                    int i=i0+a, j=j0+b; idx_t id = I(i,j);
                    double w = wx[a]*wy[b];
                    double dNx = dwx[a]*wy[b];
                    double dNy = wx[a]*dwy[b];

                    double dm = w * m;
                    double dpx = w * m * vx;
                    double dpy = w * m * vy;

                    // atomici per thread-safety
                    #pragma omp atomic
                    g.m[id] += dm;
                    #pragma omp atomic
                    g.px[id] += dpx;
                    #pragma omp atomic
                    g.py[id] += dpy;

                    // forze interne: - Vp * Pk * gradN
                    double fx_i = -(P11*dNx + P12*dNy) * Vp;
                    double fy_i = -(P21*dNx + P22*dNy) * Vp;

                    #pragma omp atomic
                    g.fx[id] += fx_i;
                    #pragma omp atomic
                    g.fy[id] += fy_i;
                }
            }
        }

        // gravità
        #pragma omp parallel for
        for(long long id=0; id<(long long)g.m.size(); ++id){
            g.fx[id] += g.m[id] * P.gx;
            g.fy[id] += g.m[id] * P.gy;
        }
    }

    void grid_update(double dt){
        // p += f dt, v = p/m; pareti rigide
        #pragma omp parallel for
        for(long long id=0; id<(long long)g.m.size(); ++id){
            if(g.m[id] > 0){
                g.px[id] += g.fx[id] * dt;
                g.py[id] += g.fy[id] * dt;
                g.vx[id] = g.px[id] / g.m[id];
                g.vy[id] = g.py[id] / g.m[id];
            } else { g.vx[id]=g.vy[id]=0.0; }
        }

        // boundary: no-penetration semplice su bordo dominio
        for(int i=0;i<P.NX;i++){
            idx_t idB = I(i,0), idT = I(i,P.NY-1);
            g.vy[idB] = std::min(0.0, g.vy[idB]);
            g.vy[idT] = std::max(0.0, g.vy[idT]);
        }
        for(int j=0;j<P.NY;j++){
            idx_t idL = I(0,j), idR = I(P.NX-1,j);
            g.vx[idL] = std::min(0.0, g.vx[idL]);
            g.vx[idR] = std::max(0.0, g.vx[idR]);
        }
    }

    void g2p_mech(double dt){
        #pragma omp parallel for
        for(long long p=0; p<(long long)(idx_t)pt->x.size(); ++p){
            double xp = pt->x[p], yp = pt->y[p];
            int i0 = (int)std::floor(xp / P.dx);
            int j0 = (int)std::floor(yp / P.dy);
            double fx = xp / P.dx - i0; if(i0<0){i0=0; fx=0;} if(i0>P.NX-2){i0=P.NX-2; fx=1;}
            double fy = yp / P.dy - j0; if(j0<0){j0=0; fy=0;} if(j0>P.NY-2){j0=P.NY-2; fy=1;}
            double wx[2] = {1.0 - fx, fx};
            double wy[2] = {1.0 - fy, fy};
            double dwx[2] = {-1.0/P.dx, 1.0/P.dx};
            double dwy[2] = {-1.0/P.dy, 1.0/P.dy};

            double vxp=0, vyp=0; // PIC
            double dvx=0, dvy=0; // per FLIP (qui non usiamo p-momentum prev, usiamo PIC/FLIP semplice)

            // gradiente velocità per aggiornare F
            double dvxx=0, dvxy=0, dvyx=0, dvyy=0;

            for(int a=0;a<2;++a){
                for(int b=0;b<2;++b){
                    int i=i0+a, j=j0+b; idx_t id = I(i,j);
                    double w = wx[a]*wy[b];
                    double dNx = dwx[a]*wy[b];
                    double dNy = wx[a]*dwy[b];
                    vxp += g.vx[id]*w; vyp += g.vy[id]*w;
                    dvxx += g.vx[id]*dNx; dvxy += g.vx[id]*dNy;
                    dvyx += g.vy[id]*dNx; dvyy += g.vy[id]*dNy;
                }
            }

            // PIC/FLIP blending: usiamo solo PIC per stabilità + pizzico di FLIP
            double vx_old = pt->dprops["vx"][p];
            double vy_old = pt->dprops["vy"][p];
            double vx_new = (1.0-P.picflip)*vx_old + P.picflip*vxp;
            double vy_new = (1.0-P.picflip)*vy_old + P.picflip*vyp;

            pt->dprops["vx"][p] = vx_new;
            pt->dprops["vy"][p] = vy_new;
            pt->x[p] += vx_new*dt;
            pt->y[p] += vy_new*dt;

            // aggiorna F: F_{n+1} = (I + dt*grad v) F_n
            Mat2 A{1 + dt*dvxx, dt*dvxy, dt*dvyx, 1 + dt*dvyy};
            Mat2 Fn = F[p];
            Mat2 Fnp1{
                A.a11*Fn.a11 + A.a12*Fn.a21,
                A.a11*Fn.a12 + A.a12*Fn.a22,
                A.a21*Fn.a11 + A.a22*Fn.a21,
                A.a21*Fn.a12 + A.a22*Fn.a22
            };
            F[p] = Fnp1;
        }
    }

    // --------- Termica ---------
    static inline double k_face(double kL,double kR){ double s=kL+kR; return s>1e-20? 2.0*kL*kR/s : 0.0; }

    void rebuild_mask_from_particles(){
        auto &mask = vars["mask"];
        std::fill(mask.begin(), mask.end(), 0.0);
        auto &phi = vars["phi"]; std::fill(phi.begin(), phi.end(), 0.0);

        #pragma omp parallel for
        for(long long p=0; p<(long long)(idx_t)pt->x.size(); ++p){
            int i = std::max(0, std::min(P.NX-1, (int)std::floor(pt->x[p]/P.dx)));
            int j = std::max(0, std::min(P.NY-1, (int)std::floor(pt->y[p]/P.dy)));
            idx_t id = I(i,j);
            #pragma omp atomic
            mask[id] += vol0[p]/(P.dx*P.dy); // frazione volumetrica
            #pragma omp atomic
            phi[id]  += pt->dprops["phi"][p];
        }
        #pragma omp parallel for
        for(long long id=0; id<(long long)mask.size(); ++id){
            mask[id] = std::clamp(mask[id], 0.0, 1.0);
            if(mask[id]>0){
                // normalizza phi media grezza dal conteggio particelle
                // (approssimazione: 1 part ~ 1 conteggio)
                phi[id] = std::clamp(phi[id], 0.0, 1.0);
            }
        }
    }

    void set_thermal_properties(){
        auto &mask = vars["mask"]; auto &phi = vars["phi"]; auto &rhoCp = vars["rhoCp"]; auto &k = vars["k"];
        #pragma omp parallel for
        for(long long id=0; id<(long long)mask.size(); ++id){
            double f = mask[id]; double ph = phi[id];
            if(f>0){
                auto h = [&](double ph){
                ph = std::clamp(ph, 0.0, 1.0);
            return ph*ph*ph*(6*ph*ph - 15*ph + 10); // 6phi^5 -15phi^4 +10phi^3
        };
             double hh = h(ph);
             double cp = (1.0-hh)*P.cp_s + hh*P.cp_l;
             double kcr= (1.0-hh)*P.k_s  + hh*P.k_l;
                rhoCp[id] = f*(P.rho*cp) + (1.0-f)*(P.rho_g*P.cp_g);
                k[id]     = f*kcr + (1.0-f)*P.k_g;
            } else {
                rhoCp[id] = P.rho_g*P.cp_g; k[id]=P.k_g;
            }
        }
    }

    double dt_mech_cfl() const{
        // stima con velocità massima e c elastica
        double vmax=1e-12; double c = std::sqrt(P.E/P.rho);
        for(idx_t p=0;p<(idx_t)pt->x.size();++p){
            double vx=pt->dprops.at("vx")[p], vy=pt->dprops.at("vy")[p];
            double v = std::hypot(vx,vy); if(v>vmax) vmax=v;
        }
        double dt1 = P.cfl_mech * std::min(P.dx,P.dy) / std::max(c, vmax+1e-9);
        return std::max(1e-8, dt1);
    }

    void thermal_step(double dt){
        auto &T=vars["T"], &phi=vars["phi"], &mask=vars["mask"], &rhoCp=vars["rhoCp"], &k=vars["k"], &S=vars["S"];
        auto &PHIDOT = vars["phidot"]; // se non l’hai aggiunto in vars, rimuovi questa riga

        auto clamp01 = [](double x){ return std::max(0.0, std::min(1.0, x)); };

        // doppio pozzo: g(phi)=phi^2(1-phi)^2  -> g'(phi)=2phi(1-phi)(1-2phi)
        auto gprime = [&](double ph){
             return 2.0*ph*(1.0-ph)*(1.0-2.0*ph);
        };

        // coupling liscio h(phi)=phi^3(6phi^2 - 15phi + 10)
        // h'(phi)=30 phi^2 (1-phi)^2
        auto hprime = [&](double ph){
           double a = ph*ph;
           double b = (1.0-ph);
           return 30.0*a*b*b;
        };

        auto idx_from = [&](int i,int j){ return (idx_t)j*P.NX + i; };

        std::vector<double> Tnew(T.size(), P.T_inf);
        std::fill(S.begin(), S.end(), 0.0);

        // -------------------- PHASE-FIELD (Allen–Cahn) + latente coerente --------------------
        std::vector<double> Tnew(T.size(), P.T_inf);
        std::vector<double> phinew(phi.size(), 0.0);
        std::vector<double> phidot(phi.size(), 0.0);
        std::fill(S.begin(), S.end(), 0.0);

        // 1) aggiorna phi con Allen–Cahn
        // dphi/dt = M * [ eps^2 Lap(phi) - W g'(phi) - lam h'(phi) (T - Tm) ]
        #pragma omp parallel for collapse(2)
        for(int j=0;j<P.NY;++j){
            for(int i=0;i<P.NX;++i){
            idx_t id = idx_from(i,j);

        // nel gas forziamo phi=0 (se vuoi farlo diffondere nel gas, togli questo if)
        if(mask[id] <= 0.0){
            phinew[id] = 0.0;
            phidot[id] = 0.0;
            continue;
        }

        double ph = clamp01(phi[id]);

        // Laplaciano con Neumann semplice (zero-gradient) via clamp sugli indici
        int il = std::max(0, i-1), ir = std::min(P.NX-1, i+1);
        int jd = std::max(0, j-1), ju = std::min(P.NY-1, j+1);

        double phL = clamp01(phi[idx_from(il,j)]);
        double phR = clamp01(phi[idx_from(ir,j)]);
        double phD = clamp01(phi[idx_from(i,jd)]);
        double phU = clamp01(phi[idx_from(i,ju)]);

        double lap =
            (phL - 2.0*ph + phR) / (P.dx*P.dx) +
            (phD - 2.0*ph + phU) / (P.dy*P.dy);

        double driveT = (T[id] - P.Tm); // >0 tende a fondere (phi->1) se lam_phi>0

        double rhs =
            (P.eps_phi*P.eps_phi) * lap
            - P.W_phi * gprime(ph)
            - P.lam_phi * hprime(ph) * driveT;

        double dphidt = P.Mphi * rhs;

        // update esplicito
        double ph_next = clamp01(ph + dt * dphidt);

        phinew[id] = ph_next;
        phidot[id] = (ph_next - ph) / dt; // per latente coerente
       }
    }

    // 2) sorgente latente: S = -rho * mask * Lm * dphi/dt
    #pragma omp parallel for
    for(long long id=0; id<(long long)S.size(); ++id){
       if(mask[id] > 0.0){
             S[id] = -(P.rho * mask[id]) * P.Lm * phidot[id];
       } else {
         S[id] = 0.0;
       }
    }

        // diffusione esplicita a k variabile (media armonica ai lati)
        auto idx_from = [&](int i,int j){ return (idx_t)j*P.NX + i; };
        #pragma omp parallel for collapse(2)
        for(int j=0;j<P.NY;++j){
            for(int i=0;i<P.NX;++i){
                idx_t id = idx_from(i,j);
                double Tij = T[id];
                double kij = k[id];
                double flux=0.0;
                if(i>0){ idx_t il=idx_from(i-1,j); flux += k_face(kij,k[il])*(T[il]-Tij)/(P.dx*P.dx); }
                if(i<P.NX-1){ idx_t ir=idx_from(i+1,j); flux += k_face(kij,k[ir])*(T[ir]-Tij)/(P.dx*P.dx); }
                if(j>0){ idx_t idn=idx_from(i,j-1); flux += k_face(kij,k[idn])*(T[idn]-Tij)/(P.dy*P.dy); }
                if(j<P.NY-1){ idx_t ids=idx_from(i,j+1); flux += k_face(kij,k[ids])*(T[ids]-Tij)/(P.dy*P.dy); }
                double rhs = flux + S[id];
                Tnew[id] = Tij + dt * rhs / std::max(1e-20, rhoCp[id]);
            }
        }

        // impone gas a T_inf se richiesto
        #pragma omp parallel for
        for(long long id=0; id<(long long)T.size(); ++id){
           T[id] = (mask[id]>0 ? Tnew[id] : (P.dirichlet_gas ? P.T_inf : Tnew[id]));
           phi[id] = phinew[id];
           PHIDOT[id] = phidot[id]; // se non esporti phidot, rimuovi
        }

        // rimappa T e phi su particelle (nearest cell)
        #pragma omp parallel for
        for(long long p=0; p<(long long)(idx_t)pt->x.size(); ++p){
            int i = std::max(0, std::min(P.NX-1, (int)std::floor(pt->x[p]/P.dx)));
            int j = std::max(0, std::min(P.NY-1, (int)std::floor(pt->y[p]/P.dy)));
            idx_t id = I(i,j);
            pt->dprops["T"][p]   = T[id];
            pt->dprops["phi"][p] = phi[id];
        }
    }

    void export_vtk(int step){
        // esporta solo variabili termiche (T, phi, mask)
        std::string fname = std::string("mpm_quadgrid_") + std::to_string(step) + ".vts";
        grid.vtk_export(fname.c_str(), vars);
    }

    void run(){
        double time=0.0;
        for(int s=0; s<P.steps; ++s){
            double dt = dt_mech_cfl();
            // --- MPM ---
            p2g_mech();
            grid_update(dt);
            g2p_mech(dt);
            // --- Termica ---
            rebuild_mask_from_particles();
            set_thermal_properties();
            thermal_step(dt);

            time += dt;
            if(s % P.output_every == 0){
                double Tavg=0, phim=0; for(idx_t p=0;p<(idx_t)pt->x.size();++p){ Tavg+=pt->dprops["T"][p]; phim+=pt->dprops["phi"][p]; }
                Tavg/=std::max<idx_t>(1,(idx_t)pt->x.size()); phim/=std::max<idx_t>(1,(idx_t)pt->x.size());
                std::cout << "step "<<s<<" t="<<time<<" dt="<<dt
                          << " T_avg="<<Tavg<<" phi_avg="<<phim<<"\n";
                export_vtk(s);
            }
        }
        std::cout << "Simulazione completata!\n";
    }
};

int main(){
    Params P;
    MPMThermoSim sim(P);
    sim.run();
    return 0;
}

