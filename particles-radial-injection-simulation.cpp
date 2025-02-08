#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <cmath>
#include <vector>
#include <random> // for generating random numbers
#include <algorithm>

//const double    v_avg = 1000;                           // [m s^-1] average velocity
//const double 	mu   = 1;                   	        // [kg s^-1 m_1] viscosity
//const double 	H    = 1.0;                   	        // [m] thickness
//const double 	KT   = .02;               	        // [J] thermal energy
//const double 	r    = 1e-2;                	        // [m] particle radius
//const double 	mob  = 1. / 6. / M_PI / mu / r; 	// [s kg^-1] mobility
//const double 	D    = KT * mob;            	        // [m^2 s^-1] diffusivity
//const int 	Np   = 2000;                    	// [-] number of particles
//const double 	dt   = 1e-2;                 	        // [s] time step
//const double 	T    = 1.0;                  	        // [s] simulation time
//const double 	g    = 9.81;                    	// [m s^-2] gravity
//const double 	m    = 0.001;                    	// [kg] particle mass
//const double    pi = 3.14;
//const double    sigma = 1.0;                            // Variance of the initial Gaussian distribution

// Function to read simulation parameters from CSV file
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

// Function to calculate fluid velocity

void fluid_velocity (double x, double y, double v_avg, double H, double& vx, double& vy){
 
 vx = 0.5 * v_avg * (1 - (y * y) / (H * H));
 vy = 0;
}

// Function to calculate particle force
void particle_force (double x, double y, double g, double m, double& fx, double& fy){

  fx = 0;
  fy = -m*g;
}

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
    
    double mob = 1.0 / (6.0 * M_PI * mu * r);
    double D = KT * mob;
    int Nt = std::ceil(T / dt);




std::vector<double> x (Np);
std::vector<double> y (Np);
std::vector<double> vx (Np);
std::vector<double> vy (Np);
std::vector<double> fx (Np);
std::vector<double> fy (Np);
std::vector<double> dxb (Np);
std::vector<double> dyb (Np);

    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
   
    // Gaussian distribution 
    std::normal_distribution<> normal(0.0, sigma); // Gaussian in x
    std::uniform_real_distribution<> uniform(0.0, 1.0); // Uniform in y
   
    // initialize particle positions
    for (int n = 0; n < Np; ++n) {  
        x[n] = normal(gen);  			// Gaussian-distributed x values
        y[n] = uniform(gen) + 0.5;              // All particles start at y = 0,5
    }
        
    // Open the CSV file for saving data
    std::ofstream file("particle_positions.csv");
    if (!file.is_open()) {
        std::cerr << "Failed to open file for writing.\n";
        return 1;
    }

    // Write the CSV header (only once)
    file << "time,x,y\n";

//Calculate the number of time steps

//int Nt = std::ceil(T/dt);
double t = 0; //initial time
int iff = 0;// File index

std::random_device rd2;  // Will be used to obtain a seed for the random number engine
std::mt19937 gen2(rd2()); // Standard mersenne_twister_engine seeded with rd()
std::normal_distribution<> disbrownian(0.0, std::sqrt (2*D*dt));
    
//Main simulation loop

  for (int it = 0; it<Nt; ++it){
  
  //Create a new CSV file for each time step
  
  // std::string particle_move = "particle_position_" + std::to_string(iff++) + ".csv";
  auto particle_move_n = std::to_string(iff++);
  std::string particle_move = "particle_position_00000.csv";
  particle_move.replace (23-particle_move_n.length (), particle_move_n.length (), particle_move_n);
  
  // Open a file to save the positions
  
  std::ofstream outFile(particle_move);
  
  // Check if the file is open
	 if (!outFile.is_open()) {
	    std::cerr << "Error opening output file." << std::endl;
	    return 1;
	 }
  
  // Write the header for the CSV file
	 outFile << "Time, Ball_X, Ball_Y \n";
  
  // Loop over all particles

    for (int n = 0; n < Np; ++n){
  
      // Compute fluid velocity and forces
  
      fluid_velocity (x[n], y[n], v_avg, H, vx[n], vy[n]);
      particle_force (x[n], y[n], g, m, fx[n], fy[n]);
  
      //Brownian motion displacements
      //std::normal_distribution<> disbrownian(0.0, std::sqrt (2*D*dt));
      
      dxb[n]=disbrownian(gen2);
      dyb[n]=disbrownian(gen2);
 
      //update particles positions
 
      x[n] = x[n] + (vx[n] + mob*fx[n]) * dt + dxb[n];
      y[n] = y[n] + (vy[n] + mob*fy[n]) * dt + dyb[n];
  
  
      // Apply boundary conditions (elastic walls)
      y[n] = std::min (H, std::max (0.0, y[n]));
      
      // Write particle position and time to the output file    
      outFile << t << "," << x[n] << "," << y[n] << "\n"; 
    
      }
      
   
   //Update time
      t += dt;
  
   //Close the time-step-specific file 
    outFile.close();
  
    }

   //Close the main output file

    file.close();
    
    std::cout << "Simulation complete. Data saved to particle_positions.csv.\n";
 
return 0;
}
