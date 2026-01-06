#include <iostream>
#include <vector>
#include <cmath>
#include <map>
#include <fstream>
#include <string>
#include <algorithm>
#include <quadgrid_cpp.h>
#include <particles.h>

using idx_t = quadgrid_t<std::vector<double>>::idx_t;

// =========================== PARAMETRI ===========================
struct Params {
    int NX = 120, NY = 120;       // celle griglia
    double Lx = 1e-3, Ly = 1e-3; // dimensioni dominio [m]
    double dt = 2e-7;           // passo temporale [s]
    int steps = 10000;           // passi totali
    int output_every = 500;     // frequenza output

    // Plasma e materiale
    double T_inf = 15000.0;     // plasma [K]
    double T0 = 300.0;          // T iniziale particella [K]
    double Tm = 2435.0;         // fusione [K]
    double rho = 5200.0;        // densità [kg/m^3]
    double cp_s = 800.0, cp_l = 1000.0;
    double k_s = 5.0, k_l = 3.0;
    double Lm = 4e5;            // calore latente [J/kg]
    double Kphi = 1e-3;         // coeff. fusione [1/s/K]
};

// ======================== SIMULATORE ============================
int main() {
    Params P;

    // --- Griglia 2D
    quadgrid_t<std::vector<double>> grid;
    grid.set_sizes(P.NY, P.NX, P.Lx / P.NX, P.Ly / P.NY);

    // --- Particelle (qui 1 sola o poche nel centro)
    std::vector<double> x, y;
    for (int j = 0; j < 5; ++j)
        for (int i = 0; i < 5; ++i)
            x.push_back(P.Lx/2 + (i-2)*P.Lx/(P.NX*2)),
            y.push_back(P.Ly/2 + (j-2)*P.Ly/(P.NY*2));

    size_t npt = x.size();
    particles_t ptcls(
        npt,
        {"id"},
        {"T","phi","mask"},
        grid, x, y
    );

    for (size_t i = 0; i < npt; ++i) {
        ptcls.dprops["T"][i] = P.T0;
        ptcls.dprops["phi"][i] = 0.0;
        ptcls.dprops["mask"][i] = 1.0; // solido
    }

    // --- Campi griglia
    std::map<std::string, std::vector<double>> vars{
        {"T", std::vector<double>(grid.num_global_nodes(), P.T_inf)},
        {"phi", std::vector<double>(grid.num_global_nodes(), 0.0)},
        {"mask", std::vector<double>(grid.num_global_nodes(), 0.0)}
    };

    // Inizializza con P2G
    ptcls.p2g(vars, {"T","phi","mask"}, {"T","phi","mask"});

    // --- Loop temporale
    for (int step = 0; step < P.steps; ++step) {
        // Diffusione semplificata (riscaldamento verso plasma)
        for (auto &T : vars["T"])
            T += P.dt * 5e3 * (P.T_inf - T); // fittizio coefficiente di scambio

        // Aggiorna fusione nelle celle solide
        for (size_t i = 0; i < vars["T"].size(); ++i) {
            double T = vars["T"][i];
            double phi = vars["phi"][i];
            if (T > P.Tm && phi < 1.0) {
                double rate = P.Kphi * (T - P.Tm) * (1.0 - phi);
                vars["phi"][i] = std::min(1.0, phi + rate * P.dt);
            }
        }

        // Aggiorna particelle
        ptcls.g2p(vars, {"T","phi"}, {"T","phi"});

        // Output periodico
        if (step % P.output_every == 0) {
            std::cout << "step " << step
                      << "  T=" << ptcls.dprops["T"][0]
                      << "  phi=" << ptcls.dprops["phi"][0] << std::endl;

            std::string fname = "melting_" + std::to_string(step) + ".vts";
            grid.vtk_export(fname.c_str(), vars);
        }
    }

    std::cout << "Simulazione completata!" << std::endl;
    return 0;
}

