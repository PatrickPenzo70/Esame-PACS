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

struct Vec2 { double x=0.0, y=0.0; };
struct MeshPoint { Vec2 position; Vec2 velocity; };

// ---------------- LETTURA CSV DELLA MESH ----------------
std::vector<MeshPoint> load_mesh_from_csv(const std::string& filename) {
    std::ifstream file(filename);
    std::vector<MeshPoint> mesh;
    if (!file.is_open()) {
        std::cerr << "Errore: impossibile aprire " << filename << "\n";
        return mesh;
    }

    std::string line;
    int x_idx=-1,y_idx=-1,vx_idx=-1,vy_idx=-1;

    if (std::getline(file, line)) {
        std::istringstream header_ss(line);
        std::string token;
        int index = 0;
        while (std::getline(header_ss, token, ',')) {
            token.erase(0, token.find_first_not_of(" \t\r\n"));
            token.erase(token.find_last_not_of(" \t\r\n") + 1);
            if (token == "Points_0") x_idx = index;
            else if (token == "Points_1") y_idx = index;
            else if (token == "Velocity_0") vx_idx = index;
            else if (token == "Velocity_1") vy_idx = index;
            index++;
        }
    }

    if (x_idx < 0 || y_idx < 0 || vx_idx < 0 || vy_idx < 0) {
        std::cerr << "Colonne CSV richieste mancanti (Points_0, Points_1, Velocity_0, Velocity_1)\n";
        return mesh;
    }

    while (std::getline(file, line)) {
        std::istringstream ss(line);
        std::vector<std::string> tokens;
        std::string val;
        while (std::getline(ss, val, ',')) tokens.push_back(val);
        if (tokens.size() <= std::max({x_idx, y_idx, vx_idx, vy_idx})) continue;

        MeshPoint mp;
        try {
            mp.position.x = std::stod(tokens[x_idx]);
            mp.position.y = std::stod(tokens[y_idx]);
            mp.velocity.x = std::stod(tokens[vx_idx]);
            mp.velocity.y = std::stod(tokens[vy_idx]);
        } catch (...) { continue; }
        mesh.push_back(mp);
    }
    return mesh;
}

// ---------------- TRASFERIMENTO MESH -> GRIGLIA (BILINEAR) ----------------
void transfer_velocity_to_grid(
    const std::vector<MeshPoint>& mesh,
    int Nx, int Ny, double dx, double dy,
    double min_x, double min_y,
    std::vector<double>& vx_grid,
    std::vector<double>& vy_grid,
    std::vector<double>& weights)
{
    if (dx <= 0.0 || dy <= 0.0) return;

    for (const auto& mp : mesh) {
        double xp = mp.position.x;
        double yp = mp.position.y;

        int i = std::clamp(int((xp - min_x) / dx), 0, Nx - 2);
        int j = std::clamp(int((yp - min_y) / dy), 0, Ny - 2);

        double sx = (xp - (min_x + i * dx)) / dx;
        double sy = (yp - (min_y + j * dy)) / dy;

        for (int di = 0; di <= 1; ++di) {
            for (int dj = 0; dj <= 1; ++dj) {
                int ni = i + di;
                int nj = j + dj;
                if (ni >= Nx || nj >= Ny) continue;

                double w = ((di == 0) ? (1 - sx) : sx) * ((dj == 0) ? (1 - sy) : sy);
                int idx = nj * Nx + ni;
                vx_grid[idx] += w * mp.velocity.x;
                vy_grid[idx] += w * mp.velocity.y;
                weights[idx] += w;
            }
        }
    }

    const double eps = 1e-12;
    for (size_t k = 0; k < vx_grid.size(); ++k) {
        if (weights[k] > eps) {
            vx_grid[k] /= weights[k];
            vy_grid[k] /= weights[k];
        } else {
            vx_grid[k] = 0.0;
            vy_grid[k] = 0.0;
        }
    }
}

// ----------------- MAIN -----------------
int main() {
    std::string filename = "Data.csv";
    auto mesh = load_mesh_from_csv(filename);
    if (mesh.empty()) {
        std::cerr << "Errore: mesh vuota o file '" << filename << "' non trovato/errato\n";
        return 1;
    }

    // Dominio
    double min_x = std::numeric_limits<double>::max();
    double max_x = std::numeric_limits<double>::lowest();
    double min_y = std::numeric_limits<double>::max();
    double max_y = std::numeric_limits<double>::lowest();
    for (const auto& mp : mesh) {
        min_x = std::min(min_x, mp.position.x);
        max_x = std::max(max_x, mp.position.x);
        min_y = std::min(min_y, mp.position.y);
        max_y = std::max(max_y, mp.position.y);
    }

    int Nx = 512, Ny = 256;
    double Lx = std::max(max_x - min_x, 1e-6);
    double Ly = std::max(max_y - min_y, 1e-6);
    double dx = Lx / double(Nx - 1);
    double dy = Ly / double(Ny - 1);

    quadgrid_t<std::vector<double>> grid;
    grid.set_sizes(Nx, Ny, dx, dy);
    size_t grid_size = grid.num_global_nodes();

    std::vector<double> vx_from_mesh(Nx * Ny, 0.0), vy_from_mesh(Nx * Ny, 0.0), w_from_mesh(Nx * Ny, 0.0);
    transfer_velocity_to_grid(mesh, Nx, Ny, dx, dy, min_x, min_y, vx_from_mesh, vy_from_mesh, w_from_mesh);

    // Trova asse del getto (massima velocità media)
    std::vector<double> vx_mean(Ny, 0.0);
    for (int j = 0; j < Ny; ++j) {
        double sum = 0.0;
        for (int i = 0; i < Nx; ++i)
            sum += vx_from_mesh[j * Nx + i];
        vx_mean[j] = sum / Nx;
    }
    int j_max = std::distance(vx_mean.begin(), std::max_element(vx_mean.begin(), vx_mean.end()));
    double y_jet_center = 0.5;
    double y_jet_width = 0.4 * (max_y - min_y);
    double y_jet_min = y_jet_center - 0.5 * y_jet_width;
    double y_jet_max = y_jet_center + 0.5 * y_jet_width;
    std::cout << "Asse del getto centrato a y = " << y_jet_center << std::endl;

    double x_jet_start = min_x;
    double x_jet_width = 0.01 * (max_x - min_x);

    // ---------------- PARTICELLE ----------------
    constexpr idx_t num_particles = 1000;
    particles_t ptcls(num_particles, {"label"}, {"m","vx","vy"}, grid);
    ptcls.dprops["m"].assign(num_particles, 1.0/double(num_particles));
    ptcls.dprops["vx"].assign(num_particles, 0.0);
    ptcls.dprops["vy"].assign(num_particles, 0.0);
    ptcls.iprops["label"].resize(num_particles);
    std::iota(ptcls.iprops["label"].begin(), ptcls.iprops["label"].end(), 0);

    // generatore casuale
    std::random_device rd;
    std::mt19937 gen(rd());

    // --- CORREZIONE: iniezione lungo l'asse centrale del getto ---
    double x_center = min_x + 0.01 * (max_x - min_x); // piccolo offset dal bordo sinistro
    double jet_radius = 0.5 * (max_y - min_y);       // spessore stretto del getto

    // Posizione media del getto
    double y_center_jet = y_jet_center;

    // Distribuzione lungo x: piccola banda in ingresso
    std::uniform_real_distribution<> dist_x(x_jet_start, x_jet_start + x_jet_width);

    // Distribuzione lungo y: concentrata attorno all’asse del getto
    std::normal_distribution<> dist_y(y_center_jet, 0.5 * (max_y - min_y)); 

    for (idx_t p = 0; p < num_particles; ++p) {
       ptcls.x[p] = dist_x(gen);
       ptcls.y[p] = std::clamp(dist_y(gen), y_jet_min, y_jet_max);
    }

    // Salva vx e vy
    {
        std::ofstream vx_out("vx_grid.csv");
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                int idx = j * Nx + i;
                vx_out << vx_from_mesh[idx];
                if (i < Nx - 1) vx_out << ",";
            }
            vx_out << "\n";
        }
    }
    {
        std::ofstream vy_out("vy_grid.csv");
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                int idx = j * Nx + i;
                vy_out << vy_from_mesh[idx];
                if (i < Nx - 1) vy_out << ",";
            }
            vy_out << "\n";
        }
    }

    // ---------------- LOOP TEMPORALE ----------------
    std::map<std::string, std::vector<double>> vars;
    vars["m"].assign(grid_size, 0.0);
    vars["vx"].assign(grid_size, 0.0);
    vars["vy"].assign(grid_size, 0.0);

    double D = 0.05;
    double dt = 1e-6;
    double sqrt2Ddt = std::sqrt(2.0 * D * dt);
    int Nt = 5000;
    int save_every = 10;
    std::normal_distribution<> brown(0.0, 1.0);

    for (int t = 0; t < Nt; ++t) {
        if (t % save_every == 0) std::cout << "Step " << t << " / " << Nt << "\n";

        std::fill(vars["m"].begin(), vars["m"].end(), 0.0);
        std::fill(vars["vx"].begin(), vars["vx"].end(), 0.0);
        std::fill(vars["vy"].begin(), vars["vy"].end(), 0.0);

        // Campo base del gas
        for (size_t ig = 0; ig < grid_size; ++ig) {
            vars["vx"][ig] += vx_from_mesh[ig];
            vars["vy"][ig] += vy_from_mesh[ig];
            vars["m"][ig] += 1e-8;
        }

        // p2g: particelle -> griglia
        ptcls.p2g(vars);

        for (size_t ig = 0; ig < grid_size; ++ig) {
            if (vars["m"][ig] > 1e-12) {
                vars["vx"][ig] /= vars["m"][ig];
                vars["vy"][ig] /= vars["m"][ig];
            }
        }

        // g2p: griglia -> particelle
        ptcls.g2p(vars);

        // Movimento particelle
        for (idx_t p = 0; p < num_particles; ++p) {
            double vx = ptcls.dprops["vx"][p];
            double vy = ptcls.dprops["vy"][p];

            ptcls.x[p] += vx * dt + brown(gen) * sqrt2Ddt;
            ptcls.y[p] += vy * dt + brown(gen) * sqrt2Ddt;

            if (ptcls.x[p] > max_x) {
                ptcls.x[p] = dist_x(gen);
                ptcls.y[p] = dist_y(gen);
                ptcls.dprops["vx"][p] = 0.0;
                ptcls.dprops["vy"][p] = 0.0;
            }

            ptcls.y[p] = std::clamp(ptcls.y[p], y_jet_min, y_jet_max);
        }

        if (t % save_every == 0) {
            double time = t * dt;
            std::ostringstream pfname;
            pfname << "particle_positions_" << std::setw(6) << std::setfill('0') << t << ".csv";
            std::ofstream pfout(pfname.str());
            pfout << "Time,label,x,y,vx,vy\n";
            for (idx_t p = 0; p < num_particles; ++p)
                pfout << time << "," << ptcls.iprops["label"][p] << ","
                      << ptcls.x[p] << "," << ptcls.y[p] << ","
                      << ptcls.dprops["vx"][p] << "," << ptcls.dprops["vy"][p] << "\n";
        }
    }

    std::cout << "✅ Simulazione completata.\n";
    return 0;
}

