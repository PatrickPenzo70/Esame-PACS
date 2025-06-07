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

using idx_t = quadgrid_t<std::vector<double>>::idx_t;

// Define a Struct for Mesh Points
struct Vec2 {
    double x = 0.0, y = 0.0;
    Vec2() = default;
    Vec2(double x_, double y_) : x(x_), y(y_) {}
};

struct MeshPoint {
    Vec2 position;
    Vec2 velocity;
};

std::vector<MeshPoint> load_mesh_from_csv(const std::string& filename) {
    std::ifstream file(filename);
    std::vector<MeshPoint> mesh;
    std::string line;

    int x_idx = -1, y_idx = -1, vx_idx = -1, vy_idx = -1;

    if (std::getline(file, line)) {
        std::istringstream header_ss(line);
        std::string token;
        int index = 0;
        while (std::getline(header_ss, token, ',')) {
            if (token == "Points_0") x_idx = index;
            else if (token == "Points_1") y_idx = index;
            else if (token == "Velocity_0") vx_idx = index;
            else if (token == "Velocity_1") vy_idx = index;
            index++;
        }
    }

    if (x_idx == -1 || y_idx == -1 || vx_idx == -1 || vy_idx == -1) {
       std::cerr << "Colonne richieste mancanti nel CSV:";
       if (x_idx == -1) std::cerr << " x";
       if (y_idx == -1) std::cerr << " y";
       if (vx_idx == -1) std::cerr << " vx";
       if (vy_idx == -1) std::cerr << " vy";
       std::cerr << '\n';
       return mesh;
       }

    while (std::getline(file, line)) {
        std::istringstream ss(line);
        std::vector<std::string> tokens;
        std::string value;

        while (std::getline(ss, value, ',')) {
            tokens.push_back(value);
        }

        if (tokens.size() <= std::max({x_idx, y_idx, vx_idx, vy_idx})) continue;

        MeshPoint mp;
        mp.position.x = std::stod(tokens[x_idx]);
        mp.position.y = std::stod(tokens[y_idx]);
        mp.velocity.x = std::stod(tokens[vx_idx]);
        mp.velocity.y = std::stod(tokens[vy_idx]);

        mesh.push_back(mp);
    }

    return mesh;
}

// Transfer Velocities to Grid
void transfer_velocity_to_grid(const std::vector<MeshPoint>& mesh,
                               int Nx, int Ny, double dx, double dy,
                               double min_x, double min_y,
                               std::vector<double>& vx_grid,
                               std::vector<double>& vy_grid,
                               std::vector<double>& weights) {
    for (const auto& mp : mesh) {
        double xp = mp.position.x;
        double yp = mp.position.y;

        // Find the lower-left corner of the grid cell
        int i = std::clamp(int((xp - min_x) / dx), 0, Nx - 2);
        int j = std::clamp(int((yp - min_y) / dy), 0, Ny - 2);

        // Local coordinates within the cell [0,1]x[0,1]
        double sx = (xp - (min_x + i * dx)) / dx;
        double sy = (yp - (min_y + j * dy)) / dy;

        // Loop over the 4 surrounding grid points
        for (int di = 0; di <= 1; ++di) {
            for (int dj = 0; dj <= 1; ++dj) {
                int ni = i + di;
                int nj = j + dj;

                if (ni >= Nx || nj >= Ny) continue;

                // Bilinear basis weight
                double w = ((di == 0) ? (1 - sx) : sx) * ((dj == 0) ? (1 - sy) : sy);

                int index = nj * Nx + ni;
                vx_grid[index] += w * mp.velocity.x;
                vy_grid[index] += w * mp.velocity.y;
                weights[index] += w;
            }
        }
    }

    // Normalize velocities by weights
    for (size_t i = 0; i < vx_grid.size(); ++i) {
        if (weights[i] > 1e-8) {
            vx_grid[i] /= weights[i];
            vy_grid[i] /= weights[i];
        }
        else {
            vx_grid[i] = 0.0;
            vy_grid[i] = 0.0;
        }
    }
}

void save_grid_velocity(const std::string& filename, int Nx, int Ny, double dx, double dy,
                        double min_x, double min_y,
                        const std::vector<double>& vx, const std::vector<double>& vy) {
    std::ofstream out(filename);
    out << "x,y,vx,vy\n";
    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            int idx = j * Nx + i;
            double x = min_x + i * dx;
            double y = min_y + j * dy;
            out << x << "," << y << "," << vx[idx] << "," << vy[idx] << "\n";
        }
    }
}

void transfer_to_grid_from_mesh(
    const std::vector<MeshPoint>& mesh,
    std::vector<double>& vx_grid,
    std::vector<double>& vy_grid,
    std::vector<double>& weight_grid,
    int Nx, int Ny,
    double dx, double dy
) 

{
    // Calculate min_x and min_y from the mesh
    double min_x = std::numeric_limits<double>::max();
    double min_y = std::numeric_limits<double>::max();
    for (const auto& mp : mesh) {
        min_x = std::min(min_x, mp.position.x);
        min_y = std::min(min_y, mp.position.y);
    }
    
    vx_grid.resize(Nx * Ny, 0.0);
    vy_grid.resize(Nx * Ny, 0.0);
    weight_grid.resize(Nx * Ny, 0.0);

    transfer_velocity_to_grid(mesh, Nx, Ny, dx, dy, min_x, min_y,
                              vx_grid, vy_grid, weight_grid);
}

int main() {
    std::string filename = "Data.csv";
    auto mesh = load_mesh_from_csv(filename);
    
    if (mesh.empty()) {
        std::cerr << "Errore: la mesh è vuota. Controlla il file '" << filename << "'\n";
        return 1;
    }
    
    // Determina bounding box
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

    int Nx = 1024;
    int Ny = 512;
    double Lx = max_x - min_x;
    double Ly = max_y - min_y;
    double dx = Lx / Nx;
    double dy = Ly / Ny;

    quadgrid_t<std::vector<double>> grid;
    grid.set_sizes(Nx, Ny, dx, dy);

    // Trasferimento velocità da mesh alla griglia
    std::vector<double> vx_from_mesh, vy_from_mesh, weight_from_mesh;
    transfer_to_grid_from_mesh(mesh, vx_from_mesh, vy_from_mesh, weight_from_mesh, Nx, Ny, dx, dy);
    
    // Pulizia: elimina eventuali valori non finiti
    for (size_t i = 0; i < vx_from_mesh.size(); ++i) {
        if (!std::isfinite(vx_from_mesh[i])) vx_from_mesh[i] = 0.0;
        if (!std::isfinite(vy_from_mesh[i])) vy_from_mesh[i] = 0.0;
    }

    save_grid_velocity("grid_velocity_initial.csv", Nx, Ny, dx, dy, min_x, min_y, vx_from_mesh, vy_from_mesh);

    constexpr idx_t num_particles = 1000;
    particles_t ptcls(num_particles, {"label"}, {"m", "vx", "vy"}, grid);
    ptcls.dprops["m"].assign(num_particles, 1.0 / static_cast<double>(num_particles));
    ptcls.dprops["vx"].resize(num_particles, 0.0);   // <-- qui la dimensione è importante!
    ptcls.dprops["vy"].resize(num_particles, 0.0);
    std::iota(ptcls.iprops["label"].begin(), ptcls.iprops["label"].end(), 0);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> normal_x(0.0, 1.0);
    std::normal_distribution<> normal_y(0.0, Ly / 3.0);
    std::normal_distribution<> brownian_noise(0.0, std::sqrt(2 * 1.0 * 1e-5)); // D=1, dt=1e-5

    // Inizializza posizioni particelle
    for (idx_t i = 0; i < num_particles; ++i) {
        ptcls.x[i] = normal_x(gen);
        ptcls.y[i] = normal_y(gen);
    }

    // Inizializza vars con dimensione corretta e valori iniziali di velocità da mesh
    size_t grid_size = grid.num_global_nodes();
    std::map<std::string, std::vector<double>> vars;
    vars["m"].resize(grid_size, 0.0);
    vars["vx"].resize(grid_size, 0.0);
    vars["vy"].resize(grid_size, 0.0);

    // Copia i dati di velocità della mesh nelle variabili della griglia
    std::copy(vx_from_mesh.begin(), vx_from_mesh.end(), vars["vx"].begin());
    std::copy(vy_from_mesh.begin(), vy_from_mesh.end(), vars["vy"].begin());

    // Debug dimensioni
    std::cout << "Dimensioni vars: m=" << vars["m"].size()
              << ", vx=" << vars["vx"].size()
              << ", vy=" << vars["vy"].size() << "\n";
    std::cout << "Numero particelle: " << num_particles << "\n";

    // Ciclo temporale
    int Nt = 2000;
    double dt = 1e-5;
    double max_vel = 1e3;

    for (int t = 0; t < Nt; ++t) {
        double time = t * dt;
        std::cout << "Time step " << t << "/" << Nt << " (t=" << time << "s)\n";

        // 1. Reset massa e velocità griglia
        std::fill(vars["m"].begin(), vars["m"].end(), 0.0);
        std::fill(vars["vx"].begin(), vars["vx"].end(), 0.0);
        std::fill(vars["vy"].begin(), vars["vy"].end(), 0.0);

        // 2. Trasferisci proprietà da particelle a griglia
        ptcls.p2g(vars);

        // 3. Normalizza velocità griglia
        for (size_t i = 0; i < vars["m"].size(); ++i) {
            if (vars["m"][i] > 1e-12) {
                vars["vx"][i] /= vars["m"][i];
                vars["vy"][i] /= vars["m"][i];
            } else {
                vars["vx"][i] = 0.0;
                vars["vy"][i] = 0.0;
            }
        }

        // 4. Trasferisci velocità da griglia a particelle
        ptcls.g2p(vars);

        // 5. Controlla valori non finiti sulle velocità particelle
        for (idx_t i = 0; i < num_particles; ++i) {
            if (!std::isfinite(ptcls.dprops["vx"][i])) ptcls.dprops["vx"][i] = 0.0;
            if (!std::isfinite(ptcls.dprops["vy"][i])) ptcls.dprops["vy"][i] = 0.0;
        }

        // 6. Debug prime particelle
        for (int i = 0; i < 5; ++i) {
            std::cout << "Particle " << i
                      << ": vx=" << ptcls.dprops["vx"][i]
                      << ", vy=" << ptcls.dprops["vy"][i] << "\n";
        }

        // 7. Aggiorna posizioni particelle con limitazioni su velocità e rumore browniano
        for (idx_t i = 0; i < num_particles; ++i) {
            double vx = std::clamp(ptcls.dprops["vx"][i], -max_vel, max_vel);
            double vy = std::clamp(ptcls.dprops["vy"][i], -max_vel, max_vel);

            ptcls.x[i] += vx * dt + brownian_noise(gen);
            ptcls.y[i] += vy * dt + brownian_noise(gen);

            ptcls.x[i] = std::clamp(ptcls.x[i], min_x, max_x);
            ptcls.y[i] = std::clamp(ptcls.y[i], min_y, max_y);

            ptcls.dprops["vx"][i] = vx;
            ptcls.dprops["vy"][i] = vy;
        }

        // 8. Salva posizioni particelle
        std::ostringstream fname;
        fname << "particle_position_" << std::setw(5) << std::setfill('0') << t << ".csv";
        std::ofstream outFile(fname.str());

        outFile << "Time,label,x,y\n";
        for (idx_t i = 0; i < num_particles; ++i) {
            outFile << time << "," << ptcls.iprops["label"][i] << "," << ptcls.x[i] << "," << ptcls.y[i] << "\n";
        }
    }

    return 0;
}



