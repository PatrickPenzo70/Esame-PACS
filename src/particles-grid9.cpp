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

using idx_t = quadgrid_t<std::vector<double>>::idx_t;

// Define a Struct for Mesh Points
struct Vec2 {
    double x, y;
};

struct MeshPoint {
    Vec2 position;
    Vec2 velocity;
};

// Load Mesh Points from CSV
std::vector<MeshPoint> load_mesh_from_csv(const std::string& filename) {
    std::ifstream file(filename);
    std::string line;
    std::vector<MeshPoint> mesh_points;

    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << "\n";
        return mesh_points;
    }

    // Skip header
    std::getline(file, line);

    while (std::getline(file, line)) {
        std::istringstream ss(line);
        std::string val;
        MeshPoint mp;

        std::getline(ss, val, ','); mp.position.x = std::stod(val);
        std::getline(ss, val, ','); mp.position.y = std::stod(val);
        std::getline(ss, val, ','); mp.velocity.x = std::stod(val);
        std::getline(ss, val, ','); mp.velocity.y = std::stod(val);

        mesh_points.push_back(mp);
    }

    return mesh_points;
}

// Transfer Velocities to Grid
void transfer_to_grid_from_mesh(
    const std::vector<MeshPoint>& mesh_points,
    std::vector<double>& vx_grid,
    std::vector<double>& vy_grid,
    std::vector<double>& weight,
    idx_t Nx, idx_t Ny, double dx, double dy
) 

  {
    for (const auto& mp : mesh_points) {
        double xp = mp.position.x;
        double yp = mp.position.y;

        int i = std::clamp(int(xp / dx), 0, static_cast<int>(Nx - 2));
        int j = std::clamp(int(yp / dy), 0, static_cast<int>(Ny - 2));

        for (int dj = 0; dj <= 1; ++dj) {
            for (int di = 0; di <= 1; ++di) {
                int ni = i + di;
                int nj = j + dj;

                if (ni < 0 || ni >= Nx || nj < 0 || nj >= Ny) continue;

                idx_t node_id = nj * Nx + ni;

                double xi = ni * dx;
                double yj = nj * dy;

                double wx = 1.0 - std::abs(xp - xi) / dx;
                double wy = 1.0 - std::abs(yp - yj) / dy;
                double w = wx * wy;

                vx_grid[node_id] += w * mp.velocity.x;
                vy_grid[node_id] += w * mp.velocity.y;
                weight[node_id] += w;
            }
        }
    }

    for (idx_t n = 0; n < vx_grid.size(); ++n) {
        if (weight[n] > 1e-12) {
            vx_grid[n] /= weight[n];
            vy_grid[n] /= weight[n];
        }
    }
}

int main(int argc, char* argv[]) {
    idx_t Nx = 512, Ny = 128;
    double Lx = 100.0;
    double dx = Lx / Nx;
    double dy = 1.0 / Ny;
    double H = 10 * Ny * dy;

    quadgrid_t<std::vector<double>> grid;
    grid.set_sizes(Nx, Ny, dx, dy);

    double v_avg = 10.0;
    double D = 1.0;
    double dt = 1e-3;
    int Nt = 10000;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> normal_x(0.0, 1.0);
    std::normal_distribution<> normal_y(0.0, H / 3.0);
    std::normal_distribution<> brownian_noise(0.0, std::sqrt(2 * D * dt));

    constexpr idx_t num_particles = 100;
    particles_t ptcls(num_particles, {"label"}, {"m", "vx", "vy"}, grid);
    ptcls.dprops["m"].assign(num_particles, 1.0 / static_cast<double>(num_particles));
    std::iota(ptcls.iprops["label"].begin(), ptcls.iprops["label"].end(), 0);

    idx_t total_nodes = grid.num_global_nodes();
    std::map<std::string, std::vector<double>> vars{
        {"m", std::vector<double>(total_nodes, 0.0)},
        {"vx", std::vector<double>(total_nodes, 0.0)},
        {"vy", std::vector<double>(total_nodes, 0.0)}
    };

    // Particle initialization
    for (idx_t i = 0; i < num_particles; ++i) {
        ptcls.x[i] = normal_x(gen);
        ptcls.y[i] = normal_y(gen);
    }

    // Load mesh and transfer velocities
    auto mesh_points = load_mesh_from_csv("Data.csv");

    std::vector<double> vx_from_mesh(total_nodes, 0.0);
    std::vector<double> vy_from_mesh(total_nodes, 0.0);
    std::vector<double> weight_from_mesh(total_nodes, 0.0);

    transfer_to_grid_from_mesh(mesh_points, vx_from_mesh, vy_from_mesh, weight_from_mesh, Nx, Ny, dx, dy);

    vars["vx"] = vx_from_mesh;
    vars["vy"] = vy_from_mesh;

    // Main time-stepping loop
    for (int t = 0; t < Nt; ++t) {
        double time = t * dt;
        std::cout << "Time step " << t << " / " << Nt << " (time = " << time << "s)\n";

        // 1. Particle to grid (mass/momentum transfer)
        ptcls.p2g(vars);

        // 2. (Optional) Apply forces or update grid velocities here

        // 3. Grid to particle (interpolate updated velocity field)
        ptcls.g2p(vars);

        // 4. Update particle positions
        for (idx_t i = 0; i < num_particles; ++i) {
            ptcls.x[i] += ptcls.dprops["vx"][i] * dt;
            ptcls.y[i] += ptcls.dprops["vy"][i] * dt;

            ptcls.x[i] = std::clamp(ptcls.x[i], 0.0, Lx);
            ptcls.y[i] = std::clamp(ptcls.y[i], -H, H);
        }

        // 5. Output particle positions
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

