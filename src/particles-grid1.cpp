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

int main (int argc, char *argv[]) {
  idx_t Nx = 64;
  idx_t Ny = 64;
  double dx = 1.0 / Nx;
  double dy = 1.0 / Ny;
  double H = Ny * dy;

  quadgrid_t<std::vector<double>> grid;
  grid.set_sizes(Nx, Ny, dx, dy);
  
  double v_avg = 1000.0;
  double D = 1.;
  double dt = 1e-3; 
  int Nt = 10000;
  int iff = 0;

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> particle_dist_x(0.0, 1.0);
  std::normal_distribution<double> particle_dist_y(0.0, H / 3.0);
  std::normal_distribution<double> brownian_noise(0.0, std::sqrt(2 * D * dt));

  std::vector<double> vx_poiseuille(grid.num_global_nodes(), 0.0);
  std::vector<double> vy_zero(grid.num_global_nodes(), 0.0);

  constexpr idx_t num_particles = 1000;
  particles_t ptcls(num_particles, {"label"}, {"m", "vx", "vy"}, grid);
  ptcls.dprops["m"].assign(num_particles, 1.0 / static_cast<double>(num_particles));
  std::iota(ptcls.iprops["label"].begin(), ptcls.iprops["label"].end(), 0);

  std::map<std::string, std::vector<double>> vars{
    {"m", std::vector<double>(grid.num_global_nodes(), 0.0)},
    {"vx", std::vector<double>(grid.num_global_nodes(), 0.0)},
    {"vy", std::vector<double>(grid.num_global_nodes(), 0.0)}
  };

  // Initial velocity field
  for (idx_t j = 0; j < Ny; ++j) {
    for (idx_t i = 0; i < Nx; ++i) {
      idx_t node_id = j * Nx + i;
      double y = j * dy;
      vx_poiseuille[node_id] = 1.5 * v_avg * (1.0 - std::pow((2.0 * y / H), 2.0)) + brownian_noise(gen);
      vy_zero[node_id] = brownian_noise(gen);
    }
  }

  // Initialize particles
  for (idx_t i = 0; i < num_particles; ++i) {
    ptcls.x[i] = std::clamp(particle_dist_x(gen), 0.0, 1.0);
    ptcls.y[i] = std::clamp(particle_dist_y(gen), -H, H);

    double y = ptcls.y[i];
    double vx_poise = 1.5 * v_avg * (1.0 - std::pow((2.0 * y / H), 2.0));

    ptcls.dprops["vx"][i] = vx_poise + brownian_noise(gen);
    ptcls.dprops["vy"][i] = brownian_noise(gen);
  }

  // Main loop
  for (int t = 0; t < Nt; ++t) {
    std::cout << "Completed timestep " << t << " / " << Nt << std::endl;

    ptcls.p2g(vars);

    std::vector<double>& vx = vars["vx"];
    std::vector<double>& vy = vars["vy"];

    // Reset velocity field each step
    for (idx_t j = 0; j < Ny; ++j) {
      for (idx_t i = 0; i < Nx; ++i) {
        idx_t id = j * Nx + i;
        double y = j * dy;
        vx[id] = 1.5 * v_avg * (1.0 - std::pow((2.0 * y / H), 2.0));
        vy[id] = 0.0;
      }
    }

    ptcls.g2p(vars);

    // Move particles
    for (idx_t i = 0; i < num_particles; ++i) {
      ptcls.x[i] += ptcls.dprops["vx"][i] * dt;
      ptcls.y[i] += ptcls.dprops["vy"][i] * dt;

      // Clamp to domain
      ptcls.x[i] = std::clamp(ptcls.x[i], 0.0, 1.0);
      ptcls.y[i] = std::clamp(ptcls.y[i], -H, H);
    }

    // Output file
    std::ostringstream fname;
    fname << "particle_position_" << std::setw(5) << std::setfill('0') << iff++ << ".csv";
    std::ofstream outFile(fname.str());
    if (!outFile.is_open()) {
      std::cerr << "Error opening output file.\n";
      return 1;
    }

    outFile << "Time,Ball_X,Ball_Y\n";
    for (idx_t i = 0; i < num_particles; ++i) {
      // Skip invalid data
      if (std::isnan(ptcls.x[i]) || std::isnan(ptcls.y[i]) ||
          std::isinf(ptcls.x[i]) || std::isinf(ptcls.y[i])) {
        continue;
      }
      outFile << t * dt << "," << ptcls.x[i] << "," << ptcls.y[i] << "\n";
    }

    outFile.close();
  }

  return 0;
}

