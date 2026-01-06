#include <fstream>
#include <functional>
#include <map>
#include <cmath>
#include <iostream>
#include <numeric>
#include <algorithm>
#include <random>
#include <particles.h>
#include <quadgrid_cpp.h>

using idx_t = quadgrid_t<std::vector<double>>::idx_t;

// -------------------- Parametri fisici --------------------
struct HeatParams {
    double k_s = 5.0, k_l = 3.0;          // W/mK
    double cp_s = 800.0, cp_l = 1000.0;   // J/kgK
    double rho = 5200.0;                  // kg/m3
    double Lm  = 4e5;                     // J/kg
    double Tm  = 2435.0;                  // K
    double T_inf = 3000.0;                // K (ambiente gas)
    double dT_enth = 100.0;               // K
    double Kphi = 1e-3;                   // 1/s/K
    double dt   = 2e-5;                   // s (CFL ok per dx=1e-5)
    double dT_clamp = 20.0;               // K per step (safety)
};

// cp apparente per entalpia
inline double cp_app(double T, const HeatParams& p) {
    if (std::fabs(T - p.Tm) <= p.dT_enth)
        return ((p.cp_s + p.cp_l)*0.5 + p.Lm/(2*p.dT_enth));
    return (T < p.Tm ? p.cp_s : p.cp_l);
}

int main() {
    // -------------------- Griglia --------------------
    const int Nx = 100, Ny = 100;
    const double hx = 1e-5, hy = 1e-5;
    quadgrid_t<std::vector<double>> grid;
    grid.set_sizes(Ny, Nx, hx, hy);

    HeatParams P;

    // -------------------- Particelle (disco) --------------------
    const idx_t num_particles = 8000;
    const double cx = 0.0005, cy = 0.0005, R = 100e-6;

    std::vector<double> x(num_particles), y(num_particles);
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> U(0.0, 1.0);
    for (idx_t i=0;i<num_particles;++i){
        double r = R*std::sqrt(U(rng)), th = 2*M_PI*U(rng);
        x[i] = cx + r*std::cos(th);
        y[i] = cy + r*std::sin(th);
    }

    // dprops: T, phi, w1 (peso=1)
    particles_t ptcls(num_particles, {"label"},
        {"T","phi","w1"}, grid, x, y);

    std::iota(ptcls.iprops["label"].begin(), ptcls.iprops["label"].end(), 0);
    ptcls.dprops["T"].assign(num_particles, 1500.0); // start caldo, ma < T_inf
    ptcls.dprops["phi"].assign(num_particles, 0.0);
    ptcls.dprops["w1"].assign(num_particles, 1.0);

    ptcls.build_mass(); // (non usata direttamente qui, ma ok averla)

    // -------------------- Variabili su griglia --------------------
    // W, Tsum, phisum servono per normalizzare P2G
    std::map<std::string, std::vector<double>> vars{
        {"W",      std::vector<double>(grid.num_global_nodes(), 0.0)},
        {"Tsum",   std::vector<double>(grid.num_global_nodes(), 0.0)},
        {"phisum", std::vector<double>(grid.num_global_nodes(), 0.0)},
        {"T",      std::vector<double>(grid.num_global_nodes(), P.T_inf)},
        {"phi",    std::vector<double>(grid.num_global_nodes(), 0.0)},
        {"k",      std::vector<double>(grid.num_global_nodes(), 0.0)},
        {"rhoCp",  std::vector<double>(grid.num_global_nodes(), 0.0)}
    };

    auto I = [&](int i,int j){ return j*Nx + i; };
    auto clampi = [&](int a,int lo,int hi){ return std::max(lo, std::min(a,hi)); };
    auto safeId = [&](int i,int j){ return I(clampi(i,0,Nx-1), clampi(j,0,Ny-1)); };

    // -------------------- Loop temporale --------------------
    double t=0.0;
    const int nsteps = 2000;

    for (int step=0; step < nsteps; ++step) {

        // (0) reset: campi di accumulo a zero; T di base = T_inf
        std::fill(vars["W"].begin(),      vars["W"].end(),      0.0);
        std::fill(vars["Tsum"].begin(),   vars["Tsum"].end(),   0.0);
        std::fill(vars["phisum"].begin(), vars["phisum"].end(), 0.0);
        std::fill(vars["T"].begin(),      vars["T"].end(),      P.T_inf);
        std::fill(vars["phi"].begin(),    vars["phi"].end(),    0.0);

        // (1) P2G: proietta pesi e quantità
        //   W   = sum N_p * w1
        //   Tsum= sum N_p * T_p
        //   phisum = sum N_p * phi_p
        ptcls.p2g(vars, std::vector<std::string>{"w1","T","phi"},
                         std::vector<std::string>{"W","Tsum","phisum"});

        // (2) normalizzazione → T, phi sui nodi con W>0; gas altrove
        for (idx_t id=0; id<vars["W"].size(); ++id){
            double W = vars["W"][id];
            if (W > 1e-12) {
                vars["T"][id]   = vars["Tsum"][id]   / W;
                vars["phi"][id] = std::clamp(vars["phisum"][id] / W, 0.0, 1.0);
            } else {
                vars["T"][id]   = P.T_inf;
                vars["phi"][id] = 0.0;
            }
        }

        // (3) proprietà apparenti
        for (idx_t id=0; id<vars["T"].size(); ++id){
            double T   = vars["T"][id];
            double phi = vars["phi"][id];
            double cp  = cp_app(T, P);
            vars["rhoCp"][id] = P.rho * cp;
            vars["k"][id]     = (1.0-phi)*P.k_s + phi*P.k_l;
        }

        // (4) diffusione esplicita + bordi Dirichlet
        std::vector<double> Tnew = vars["T"];
        const auto& k = vars["k"];
        const auto& rhoCp = vars["rhoCp"];

        double max_dT = 0.0;

        for (int j=0;j<Ny;++j){
            for (int i=0;i<Nx;++i){
                int id = I(i,j);

                // Bordo: Dirichlet a T_inf
                if (i==0 || j==0 || i==Nx-1 || j==Ny-1) {
                    Tnew[id] = P.T_inf;
                    continue;
                }

                double T = vars["T"][id];

                // 5-points con k variabile (semplificato; media aritmetica)
                int idL = I(i-1,j), idR = I(i+1,j), idD = I(i,j-1), idU = I(i,j+1);
                double flux =
                    k[idL]*(vars["T"][idL] - T) +
                    k[idR]*(vars["T"][idR] - T) +
                    k[idD]*(vars["T"][idD] - T) +
                    k[idU]*(vars["T"][idU] - T);

                flux /= (hx*hx);

                double dT = P.dt * flux / (rhoCp[id] + 1e-12);
                dT = std::clamp(dT, -P.dT_clamp, P.dT_clamp); // safety
                max_dT = std::max(max_dT, std::fabs(dT));

                double Tn = T + dT;
                if (!std::isfinite(Tn)) Tn = P.T_inf;
                Tnew[id] = std::clamp(Tn, 0.0, 2.0*P.T_inf);
            }
        }

        // (5) aggiorna T
        vars["T"].swap(Tnew);

        // (6) aggiorna φ SOLO sui nodi con particelle (W>0)
        for (idx_t id=0; id<vars["phi"].size(); ++id){
            if (vars["W"][id] > 1e-12 && vars["T"][id] > P.Tm && vars["phi"][id] < 1.0){
                double dphi = P.Kphi * (vars["T"][id] - P.Tm) * (1.0 - vars["phi"][id]) * P.dt;
                dphi = std::min(dphi, 0.01); // safety
                vars["phi"][id] = std::clamp(vars["phi"][id] + dphi, 0.0, 1.0);
            }
        }

        // (7) G2P: rimappa T e φ alle particelle
        ptcls.g2p(vars, std::vector<std::string>{"T","phi"},
                        std::vector<std::string>{"T","phi"});

        // (8) diagnostica + output
        if (step % 100 == 0) {
            double Tavg=0.0, phim=0.0;
            for (idx_t p=0;p<num_particles;++p){ Tavg+=ptcls.dprops["T"][p]; phim+=ptcls.dprops["phi"][p]; }
            Tavg/=num_particles; phim/=num_particles;

            auto mmT   = std::minmax_element(vars["T"].begin(), vars["T"].end());
            auto mmphi = std::minmax_element(vars["phi"].begin(), vars["phi"].end());
            std::cout << "step " << step
                      << " t=" << t
                      << " Tavg=" << Tavg
                      << " Tmin=" << *mmT.first
                      << " Tmax=" << *mmT.second
                      << " phi_avg=" << phim
                      << " phi_max=" << *mmphi.second
                      << " max_dT=" << max_dT
                      << std::endl;

            // CSV particelle
            std::ofstream f("heat_out_"+std::to_string(step)+".csv");
            f << "x,y,T,phi\n";
            for (idx_t p=0;p<num_particles;++p)
                f << ptcls.x[p] << "," << ptcls.y[p] << ","
                  << ptcls.dprops["T"][p] << "," << ptcls.dprops["phi"][p] << "\n";
            f.close();
        }

        t += P.dt;
    }

    return 0;
}

