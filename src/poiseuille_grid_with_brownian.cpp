#include <algorithm>
#include <random>
#include <quadgrid_cpp.h>
#include <particles.h>
#include <map>
#include <iostream>

using idx_t = quadgrid_t<std::vector<double>>::idx_t;

int main(int argc, char *argv[]) {
    // Create the grid
    quadgrid_t<std::vector<double>> grid;
    grid.set_sizes(32, 32, 1./32., 1./32.);
    
    constexpr idx_t num_particles = 1000000;
    particles_t ptcls(num_particles, {"label"}, {"m", "vx", "vy"}, grid);
    
    // Set uniform mass for all particles
    ptcls.dprops["m"].assign(num_particles, 1. / static_cast<double>(num_particles));
    
    // Assign unique labels to particles
    idx_t ilabel = 0;
    std::iota(ptcls.iprops["label"].begin(), ptcls.iprops["label"].end(), ilabel);
    
    // Random number generator for Brownian noise
    std::random_device rd;
    std::mt19937 gen(rd());  // Mersenne Twister RNG
    std::normal_distribution<double> brownian_noise(0.0, 0.2);  // Mean=0, StdDev=0.2
    
    // Poiseuille flow parameters
    double Vmax = 1500.0;  // Maximum velocity at the center
    double h = 0.5;     // Half of channel height
    
    // Assign initial velocities with Poiseuille profile + Brownian noise
    for (idx_t i = 0; i < num_particles; ++i) {
        double y = static_cast<double>(rand()) / RAND_MAX;  // Random y in [0, 1]
        double vx = Vmax * (1 - (y / h) * (y / h)) + brownian_noise(gen);
        double vy = brownian_noise(gen);  // Brownian perturbation in y

        ptcls.dprops["vx"][i] = vx;
        ptcls.dprops["vy"][i] = vy;
    }
    
    // Create storage for grid variables
    std::map<std::string, std::vector<double>> vars{
        {"m", std::vector<double>(grid.num_global_nodes(), 0.0)},
        {"vx", std::vector<double>(grid.num_global_nodes(), 0.0)},
        {"vy", std::vector<double>(grid.num_global_nodes(), 0.0)}
    };
    
    // Transfer particle properties to the grid
    ptcls.p2g(vars);

    // Print mass at each grid node
    for (auto ii : vars["m"]) {
        std::cout << ii << std::endl;
    }

    return 0;
}
    
    
    
    
    
    
