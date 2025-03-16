#include <iostream>        // For input and output operations
#include <fstream>         // For file handling (reading and writing files)
#include <sstream>         // For string stream operations
#include <unordered_map>   // For using hash maps (storing parameters)
#include <cmath>           // For mathematical functions like sqrt, clamp, etc.
#include <vector>          // For dynamic arrays (vectors)
#include <random>          // For random number generation
#include <algorithm>       // For various algorithms (like std::for_each)
#include <exception>       // For handling exceptions
#include <execution>       // For parallel execution (in std::for_each)
#include <functional>      // For using std::function (function pointers or lambdas)
#include "counter.h"       // This is a custom header file for range counting



// This function reads a CSV file and stores parameters in a map (unordered_map):
// File handling: It opens the CSV file, reads each line, and splits it into key-value pairs (parameter name and value).
// Error handling: If the file cannot be opened, it prints an error message and returns false.

bool readParameters(const std::string& filename, std::unordered_map<std::string, double>& params) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << "\n";
        return false;
    }

    std::string line, key;
    double value;

// Skip header
    std::getline(file, line);

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        if (std::getline(ss, key, ',') && ss >> value) {
            params[key] = value;
        }
    }

    file.close();
    return true;
}

// QuadGrid class to efficiently manage particle positions in a Cartesian grid
// QuadGrid class:
// The QuadGrid class is used to manage particle positions in a 2D grid:
// Grid representation: It divides the simulation space into a grid and stores particle IDs in each grid cell.
// Insertion method: Particles are inserted into a cell based on their position.
// Clear method: Clears the grid at the start of each time step to avoid mixing old and new positions.

class QuadGrid {
private:
    int grid_size_x, grid_size_y;
    double cell_width, cell_height;
    std::vector<std::vector<std::vector<int>>> grid;

public:
    QuadGrid(int size_x, int size_y, double width, double height)
        : grid_size_x(size_x), grid_size_y(size_y), cell_width(width), cell_height(height),
          grid(size_x, std::vector<std::vector<int>>(size_y)) {}

    void insert(int particle_id, double x, double y) {
        int i = std::min(grid_size_x - 1, std::max(0, static_cast<int>(x / cell_width)));
        int j = std::min(grid_size_y - 1, std::max(0, static_cast<int>(y / cell_height)));
        grid[i][j].push_back(particle_id);
    }

    void clear() {
        for (auto &col : grid) {
            for (auto &cell : col) {
                cell.clear();
            }
        }
    }
};


// Function to calculate fluid velocity at a given position (x, y): it computes the velocity vx based on the position of the particle relative to the height H and the average fluid velocity v_avg.

void fluid_velocity (double x, double y, double v_avg, double mu, double H, double& vx, double& vy){
 
 vx = 1.5 * v_avg * (1 - y*y/(H*H)); // Velocity decreases with y^2
 vy = 0; // There is no velocity component in the y-direction
}

// particle_force function:
// This function calculates the force on a particle due to gravity:
// The gravitational force fy is calculated using the particle's mass m and gravity g.

void particle_force (double x, double y, double g, double m, double& fx, double& fy){

  fx = 0; // No force in x direction
  fy = -m*g; // Gravitational force acting downward
}

// stepper class responsible for updating the positions of the particles using the equations of motion:
// Input parameters: The class takes parameters such as particle positions (x, y), fluid velocity, forces, and Brownian motion characteristics (through normal()).
// Update: For each particle, the class updates its position based on the forces and velocity (including Brownian motion).
// Boundary conditions: The class also handles boundary conditions (elastic walls), ensuring particles don't move beyond the simulation boundaries.

class stepper {

private:
  std::vector<double> &x;
  std::vector<double> &y;
  double v_avg, mu, H, KT, r, mob, D, dt, g, m; 
  std::function<double (void)> normal;
  
public:
  
  stepper (std::vector<double> &x_, std::vector<double> &y_, double v_avg_, double mu_,
	   double H_, double KT_, double r_, double mob_, double D_, double dt_, double g_, 
	   double m_, std::function<double (void)> &noise_)
    : x(x_), y(y_), v_avg{v_avg_}, mu{mu_}, H{H_}, KT{KT_}, r{r_}, mob{mob_}, D{D_}, dt{dt_}, g{g_},
      m{m_}, normal(noise_) { };
  
  void operator() (int n) {
    double vx, vy, fx, fy, dxb, dyb;
 
    // Compute fluid velocity and forces
  
    fluid_velocity (x[n], y[n], v_avg, mu, H, vx, vy);
    particle_force (x[n], y[n], g, m, fx, fy);
  
    //Brownian motion displacements
    
    dxb = normal() * std::sqrt (2*D*dt);
    dyb = normal() * std::sqrt (2*D*dt);
 
    //update particles positions
    x[n] += (vx + mob*fx) * dt + dxb;
    y[n] += (vy + mob*fy) * dt + dyb;
  
     // Apply boundary conditions (elastic walls)
    y[n] = std::clamp(y[n], -H, H);
   }; 

   void save(const std::string &filename) const {
     std::ofstream outFile (filename);
     int Np = x.size ();

    // Check if the file is open
	 if (!outFile.is_open()) {
	    std::cerr << "Error opening output file." << std::endl;
	 }
  
    // Write the header for the CSV file
	 outFile << "Ball_X, Ball_Y \n";
	 
    // Write particle position to the output file 
    for (int n = 0; n < Np; ++n){   
      outFile << x[n] << "," << y[n] << "\n"; 
      }

    //Close the time-step-specific file 
    outFile.close();
   }
  };


// Main Function:
// Parameter reading: Reads simulation parameters from a CSV file using readParameters.
// Random initialization: Initializes particle positions randomly within the grid.
// Simulation loop: Performs the simulation over Nt time steps:
// Particles' positions are updated using the stepper class.
// The QuadGrid class to manage particle positions within a grid, which could help with spatial optimization if further code (like interaction calculations) was added.
// Particle positions are saved every 1000 steps.
// Parallel execution using std::for_each with the std::execution::par policy to parallelize updates for each particle.

int main() {

    std::unordered_map<std::string, double> params;

    // Read parameters from CSV file
    if (!readParameters("Data.csv", params)) {
        return 1;
    }

    // Assign values from the map

    double v_avg = params["v_avg"];
    int Np = static_cast<int>(params["Np"]);
    double dt = params["dt"];
    double T = params["T"];
    double g = params["g"];
    double m = params["m"];
    double H = params["H"];
    double KT = params["KT"];
    double r = params["r"];
    double sigma = params["sigma"];
    double mu = params["mu"];
    double mob = 1.0;
    double D = KT * mob;
    int Nt = std::ceil(T / dt);
    int size_x = params["size_x"];
    int size_y = params["size_y"];
    double width = params["width"];
    double height = params["height"];
    

    std::vector<double> x (Np);
    std::vector<double> y (Np);
    
    
    // Initialize particle positions randomly
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::normal_distribution<> normal_x(0.0, 1.0); // Mean = 0, Std Dev = 1
    std::uniform_real_distribution<> uniform_y(0.0, 0.125); // y in [0.5, 0.625]
    
   
    for (int n = 0; n < Np; ++n) {  
        x[n] = normal_x(gen);  	        // Gaussian-distributed x values
        y[n] = uniform_y(gen);            // All particles start at y = 0,5
    }
        
    // Brownian motion noise function
        
    std::random_device rd2;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen2(rd2()); // Standard mersenne_twister_engine seeded with rd() 
    
    std::function<double ()> noise = [&gen2, &normal_x] (){
         return normal_x(gen2);
    };

    stepper state (x, y, v_avg, mu, H, KT, r, mob, D, dt, g, m, noise);
    range counter (0, Np);
    QuadGrid grid(size_x, size_y, width, height);
    
    
      // Main simulation loop
    
    for (int it = 0; it < Nt; ++it) {
        grid.clear();
        for (int n = 0; n < Np; ++n) {
            state(n);
            grid.insert(n, x[n], y[n]);
        }
    }

    //Calculate the number of time steps
    double t = 0; //initial time
    int iff = 0;// File index
    
    //Main simulation loop
    bool printframe = false;
    std::ofstream outFile;
    for (int it = 0; it<Nt; ++it){

    // Loop over all particles
    std::for_each (std::execution::par, counter.begin (), counter.end (), state);
    
    printframe = false;
    if ((it % 20) == 0) printframe = true;
    
    
    if (printframe) {
     //Create a new CSV file for each time step
     auto particle_move_n = std::to_string(iff++);
     std::string particle_move = "particle_position_00000.csv";
     particle_move.replace (23-particle_move_n.length (), particle_move_n.length (), particle_move_n); state.save (particle_move);
    }
    
       
   //Update time
      t += dt;
    }
    
    std::cout << "Simulation complete. Data saved to particle_positions.csv.\n";
    return 0;
   }
   
   
   
   
