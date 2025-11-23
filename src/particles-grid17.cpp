#include <algorithm>
#include <random>
#include <quadgrid_cpp.h>
#include <particles.h>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <cmath>
#include <vector>
#include <iomanip>
#include <limits>
#include <string>
#include <numeric>

using idx_t = quadgrid_t<std::vector<double>>::idx_t;

// ---------------- SALVATORE CSV ----------------
void save_particles_csv(const particles_t &ptcls, double time, int it) {
    std::ostringstream pfname;
    pfname << "particle_positions_" << std::setw(6) << std::setfill('0') << it << ".csv";
    std::ofstream pfout(pfname.str());
    pfout << "Time,label,x,y,vx,vy\n";
    for (idx_t p = 0; p < ptcls.num_particles; ++p)
        pfout << time << "," << ptcls.iprops.at("label")[p] << ","
              << ptcls.x[p] << "," << ptcls.y[p] << ","
              << ptcls.dprops.at("vx")[p] << "," << ptcls.dprops.at("vy")[p] << "\n";
}

// ---------------- PARTICLE TO GRID (P2G) ----------------
void p2g(const particles_t &ptcls,
         std::vector<double> &mass_grid,
         std::vector<double> &momx_grid,
         std::vector<double> &momy_grid,
         int Nx, int Ny, double dx, double dy)
{
    std::fill(mass_grid.begin(), mass_grid.end(), 0.0);
    std::fill(momx_grid.begin(), momx_grid.end(), 0.0);
    std::fill(momy_grid.begin(), momy_grid.end(), 0.0);

    for (idx_t p = 0; p < ptcls.num_particles; ++p) {
        double xp = ptcls.x[p];
        double yp = ptcls.y[p];
        double mp = ptcls.dprops.at("m")[p];
        double vxp = ptcls.dprops.at("vx")[p];
        double vyp = ptcls.dprops.at("vy")[p];

        int i = std::clamp(int(xp / dx), 0, Nx - 2);
        int j = std::clamp(int(yp / dy), 0, Ny - 2);
        double fx = (xp - i * dx) / dx;
        double fy = (yp - j * dy) / dy;

        double w[2][2] = {
            {(1 - fx) * (1 - fy), (1 - fx) * fy},
            {fx * (1 - fy), fx * fy}
        };

        for (int di = 0; di < 2; ++di) {
            for (int dj = 0; dj < 2; ++dj) {
                int gi = i + di;
                int gj = j + dj;
                if (gi >= 0 && gi < Nx && gj >= 0 && gj < Ny) {
                    int idx = gj * Nx + gi;
                    double weight = w[di][dj];
                    mass_grid[idx] += mp * weight;
                    momx_grid[idx] += mp * vxp * weight;
                    momy_grid[idx] += mp * vyp * weight;
                }
            }
        }
    }
}

// ---------------- GRID TO PARTICLE (G2P) ----------------
void g2p(particles_t &ptcls,
         const std::vector<double> &vx_grid,
         const std::vector<double> &vy_grid,
         int Nx, int Ny, double dx, double dy)
{
    for (idx_t p = 0; p < ptcls.num_particles; ++p) {
        double xp = ptcls.x[p];
        double yp = ptcls.y[p];
        int i = std::clamp(int(xp / dx), 0, Nx - 2);
        int j = std::clamp(int(yp / dy), 0, Ny - 2);

        double fx = (xp - i * dx) / dx;
        double fy = (yp - j * dy) / dy;

        double w[2][2] = {
            {(1 - fx) * (1 - fy), (1 - fx) * fy},
            {fx * (1 - fy), fx * fy}
        };

        double vx_interp = 0.0, vy_interp = 0.0;

        for (int di = 0; di < 2; ++di) {
            for (int dj = 0; dj < 2; ++dj) {
                int gi = i + di;
                int gj = j + dj;
                if (gi >= 0 && gi < Nx && gj >= 0 && gj < Ny) {
                    int idx = gj * Nx + gi;
                    double weight = w[di][dj];
                    vx_interp += vx_grid[idx] * weight;
                    vy_interp += vy_grid[idx] * weight;
                }
            }
        }

        ptcls.dprops["vx"][p] = vx_interp;
        ptcls.dprops["vy"][p] = vy_interp;
            }
        }

// ---------------- MAIN ----------------
int main() {
    // Dominio
    double min_x = 0.0, max_x = 1.0;
    double min_y = 0.0, max_y = 1.0;
    int Nx = 256, Ny = 128;
    double dx = (max_x - min_x) / (Nx - 1);
    double dy = (max_y - min_y) / (Ny - 1);

    quadgrid_t<std::vector<double>> grid;
    grid.set_sizes(Nx, Ny, dx, dy);
    size_t grid_size = grid.num_global_nodes();

    // --- CAMPO DEL GETTO centrato a y=0.5 ---
    std::vector<double> vx_from_mesh(Nx * Ny, 0.0), vy_from_mesh(Nx * Ny, 0.0);
    double y_jet_center = 0.5;
    for (int j = 0; j < Ny; ++j) {
        double y = min_y + j * dy;
        double weight = std::exp(-std::pow((y - y_jet_center) / 0.1, 2));
        for (int i = 0; i < Nx; ++i) {
            int idx = j * Nx + i;
            vx_from_mesh[idx] = 20.0 * weight;
            vy_from_mesh[idx] = 0.0;
        }
    }

    // --- PARTICELLE ---
    constexpr idx_t num_particles = 10000;
    particles_t ptcls(num_particles, {"label"}, {"m","vx","vy","ax","ay"}, grid);
    ptcls.dprops["m"].assign(num_particles, 1e1 / double(num_particles));
    ptcls.iprops["label"].resize(num_particles);
    std::iota(ptcls.iprops["label"].begin(), ptcls.iprops["label"].end(), 0);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist_x(0.0, 0.05);
    std::normal_distribution<> dist_y(y_jet_center, 0.02);

    for (idx_t p = 0; p < num_particles; ++p) {
        ptcls.x[p] = dist_x(gen);
        ptcls.y[p] = std::clamp(dist_y(gen), 0.0, 1.0);
    }

    for (idx_t p = 0; p < num_particles; ++p) {
        int i = std::clamp(int(ptcls.x[p] / dx), 0, Nx - 2);
        int j = std::clamp(int(ptcls.y[p] / dy), 0, Ny - 2);
        int idx = j * Nx + i;
        ptcls.dprops["vx"][p] = vx_from_mesh[idx];
        ptcls.dprops["vy"][p] = 0.0;
    }

    // --- PARAMETRI FISICI ---
    double D = 0.02;
    double dt = 1.0e-3;
    // particle physical size and material properties
    double d = 240e-2;     // particle diameter [m]
    double rho_p = 5220;  // particle density [kg/m^3]
    double mu = 1.8e-5;   // dynamic viscosity of air [Pa·s]

    // Stokes relaxation time for spheres
    double tau = (rho_p * d * d) / (18.0 * mu);
    double sqrt2Ddt = std::sqrt(2.0 * D * dt);
    int Nt = 20000;
    int save_every = 10;
    std::normal_distribution<> brown(0.0, 1.0);

    std::vector<double> vgas_x = vx_from_mesh;
    std::vector<double> vgas_y = vy_from_mesh;
    std::vector<double> mass_grid(Nx * Ny, 0.0);
    std::vector<double> momx_grid(Nx * Ny, 0.0);
    std::vector<double> momy_grid(Nx * Ny, 0.0);

    // --- LOOP TEMPORALE ---
    for (int it = 0; it < Nt; ++it) {
        // --- P2G ---
        p2g(ptcls, mass_grid, momx_grid, momy_grid, Nx, Ny, dx, dy);

        // Calcolo velocità media sulla griglia
        for (size_t i = 0; i < mass_grid.size(); ++i) {
            if (mass_grid[i] > 0.0) {
                vgas_x[i] = momx_grid[i] / mass_grid[i];
                vgas_y[i] = momy_grid[i] / mass_grid[i];
            } else {
                vgas_x[i] = 0.0;
                vgas_y[i] = 0.0;
            }
        }
        
        // store old particle velocity (needed for relaxation)
        auto vx_old_vec = ptcls.dprops["vx"];
        auto vy_old_vec = ptcls.dprops["vy"];
        
        // --- G2P ---
        g2p(ptcls, vgas_x, vgas_y, Nx, Ny, dx, dy);
        
        // --- Relaxation (Stokes drag) ---
           double r = std::exp(-dt / tau);          // relaxation factor
           auto &vx = ptcls.dprops["vx"];
           auto &vy = ptcls.dprops["vy"];

        for (idx_t p = 0; p < num_particles; ++p) {
           double ux = vx[p];                // fluid/grid velocity after G2P
           double uy = vy[p];

           double vx_old = vx_old_vec[p];    // particle velocity BEFORE G2P
           double vy_old = vy_old_vec[p];

           vx[p] = ux + (vx_old - ux) * r;   // exact dv/dt = -(v - u)/tau
           vy[p] = uy + (vy_old - uy) * r;
        }
        
        // Gravity + buoyancy
        double g = 9.81;
        double rho_f = 1.6;       // gas density aria/argon approx.
        double coeff = (rho_p - rho_f) / rho_p;
        for (idx_t p = 0; p < num_particles; ++p) {
             vy[p] -= g * coeff * dt;
        }
 
        // --- Dinamica particellare ---
        for (idx_t p = 0; p < num_particles; ++p) {
             double rx = brown(gen) * sqrt2Ddt;
             double ry = brown(gen) * sqrt2Ddt;

             ptcls.x[p] += ptcls.dprops["vx"][p] * dt + rx;
             ptcls.y[p] += ptcls.dprops["vy"][p] * dt + ry;

       // ---- Reinjection with memory loss ----
       if (ptcls.x[p] < 0.0) {
           ptcls.x[p] = 0.0;    // do NOT reset velocity here
       }

       if (ptcls.x[p] > max_x) {
           ptcls.x[p] = 0.0;
           ptcls.y[p] = std::clamp(dist_y(gen), 0.0, 1.0);

           int j = std::clamp(int(ptcls.y[p] / dy), 0, Ny - 2);
           int idx = j * Nx + 0;

           ptcls.dprops["vx"][p] = vx_from_mesh[idx];
           ptcls.dprops["vy"][p] = vy_from_mesh[idx];
       }

        ptcls.y[p] = std::clamp(ptcls.y[p], 0.0, 1.0);
    }
        // Output periodico
        if (it % save_every == 0) {
            double time = it * dt;
            std::cout << "Step " << it << "/" << Nt << std::endl;
            save_particles_csv(ptcls, time, it);
        }
    }

    std::cout << "\n✅ Simulazione completata. File CSV generati.\n";
    return 0;
}

