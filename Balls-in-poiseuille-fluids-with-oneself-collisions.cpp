// MOTION SIMULATING OF BALLS/PARTICLES IN A POISEUILLE FLUID FLOW IN A PIPELINE.
// POISEUILLE FLUID FLOW GENERALLY REFERS TO THE FLOW OF A VISCOUS FLUID
// GOVERNED BY THE NAVIER-STOKES EQUATIONS FOR LAMINAR FLOW.
// THE VELOCITY PROFILE OF A POISEULLE FLOW IS TYPICALLY PARABILIC, WITH MAXIMUM VELOCITY 
// AT THE CENTER OF THE PIPE AND ZERO AT THE BOUNDARIES DUE TO NO-SLIP CONDITIONS.

// THE FLUID FLOW WILL AFFECT THE BALLS, CAUSING THEM TO MOVE ACCORDING TO DRAG FROCES AND GRAVITY 
// 
// IT COMPUTES THEIR POSITIONS OVER THE TIME.


// 1) POISEUILLE FLOW DEFINITION: The velocity of the fluid at a point (x.y) is parabolic in the 
//    cross-section direction (let's say along the y-axis for a pipe flow). The velocity u(y) could be defined as:
//    u(y) = Umax(1 - y*y/R*R)

// Where:
// Umax is the maximum velocity at the center of the pipe,
// 
// R is the radius of the pipe 
// 
// y is the position along the cross-section direction of the pipe
//
// 2) BALL'S MOTION: A particle in the fluid will experience a drag force based on the local fluid 
//    velocity and the particle's velocity. This can be modeled using Stokes' drag law for 
//    low Reynolds numbers:

//    F_drag_force = 6_M_PI * mu * r (v_fluid - v_particle)

// Where:
// mu is the dynamic viscosity of the fluid

// r is the radius of the particle

// v_fluid is the velocity of the fluid at the particle's position

// v_particle is the velocity of the particle

// 3) NEWTON'S SECOND LAW: The particle's motion can be computed using Newton's second law:

// m*dv_particle / dt = F_drag_force + F_gravity

// Where: 
// m is the mass of the particle,

// F_gravity = m * g is the gravitational force acting on the particle


// 4) NUMERICAL INTEGRATION: To simulate the motion of the particle over the time, I can use a 
//    time-stepping method such as Euler or Verlet integration to update the position and velocity
//    of the particle at each time step.

// 5) To save the position of the particles

// 6) Collision Detection: We check whether two particles overlap, i.e., 
//    if the distance between their centers is less than or equal to the sum of their radii.


#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <random> // for generating random numbers
#include <iomanip>
#include <algorithm>



/*
G   = 30;                  %% [kg s^-2 m^-2] head
mu  = 1;                   %% [kg s^-1 m_1] viscosity
h   = 1;                   %% [m] thickness
KT  = .0001;                 %% [J] thermal energy
r   = 1e-2;                %% [m] particle radius
mob = 1 / 6 / pi / mu / r; %% [s kg^-1] mobility
D   = KT * mob;            %% [m^2 s^-1] diffusivity
Np = 2000;                 %% [-] number of particles
dt = 1e-3;                 %% [s] time step
T  = 300;                  %% [s] simulation tim
g  = 0;                    %% [m s^-2] gravity
m  = 1;                    %% [kg] particle mass
*/

/*%% fluid velocity profile
function [vx, vy] = fluid_velocity (x, y, G, mu, h)
  %% Stationary Couette or
  %% plane Poiseuille flow
  vx = G/mu/2 .* y .* (h - y);
  vy = 0;
endfunction */

// Function to calculate fluid velocity

void fluid_velocity (double x, double y, double G, double mu, double h, double& vx, double& vy){
 
 vx = G/mu/2 * y * (h - y);
 vy = 0;

}




/*
function [fx, fy] = particle_force (x, y, g, m)
  fx = 0;
  fy = -m*g;
endfunction

*/


// Function to calculate particle force

void particle_force (double x, double y, double g, double m, double& fx, double& fy){

  fx = 0;
  fy = -m*g;
  
}





int main() {

double 	G   = 30;                  	// [kg s^-2 m^-2] head
double 	mu  = 1;                   	// [kg s^-1 m_1] viscosity
double 	h   = 1;                   	// [m] thickness
double 	KT  = .0001;               	// [J] thermal energy
double 	r   = 1e-2;                	// [m] particle radius
double 	mob = 1. / 6. / M_PI / mu / r; 	// [s kg^-1] mobility
double 	D   = KT * mob;            	// [m^2 s^-1] diffusivity
int 	Np = 2000;                    	// [-] number of particles
double 	dt = 1e-2;                 	// [s] time step
double 	T  = 1.0;                  	// [s] simulation time
double 	g  = 0;                    	// [m s^-2] gravity
double 	m  = 1;                    	// [kg] particle mass

/*
x  = randn (Np, 1) * 10;
y  = randn (Np, 1)/8 + 0.5;
*/

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
    std::uniform_real_distribution<> disx(0.0, 0.1);
    std::uniform_real_distribution<> disy(0.0, 0.1);
    
    // initialize particle positions
    for (int n = 0; n < Np; ++n) {  
        x[n]=disx(gen);
        y[n]=disy(gen);
    }
        

// Open the CSV file for saving data
    std::ofstream file("particle_positions.csv");
    if (!file.is_open()) {
        std::cerr << "Failed to open file for writing.\n";
        return 1;
    }


   // Write the CSV header (only once)
    file << "time,x,y\n";



/*
Nt = ceil (T/dt);
t  = 0;
iff = 0;
*/

//Calculate the number of time steps

int Nt = std::ceil(T/dt);

double t = 0; //initial time

int iff = 0;// File index

/*
for it = 1 : Nt

  [vx, vy] = fluid_velocity (x, y, G, mu, h);
  [fx, fy] = particle_force (x, y, g, m);
  
  dxb = randn (Np, 2) * sqrt (2*D*dt);
  x = x + (vx + mob*fx) * dt + dxb(:, 1);
  y = y + (vy + mob*fy) * dt + dxb(:, 2);

  %% upper and lower (unelastic) wall conditions
  y = min (h, max (0, y));
  t = t + dt;

*/


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
  
      fluid_velocity (x[n], y[n], G, mu, h, vx[n], vy[n]);
      particle_force (x[n], y[n], g, m, fx[n], fy[n]);
  
      //Brownian motion displacements
      //std::normal_distribution<> disbrownian(0.0, std::sqrt (2*D*dt));
      
      dxb[n]=disbrownian(gen2);
      dyb[n]=disbrownian(gen2);
 
      //update particles positions
 
      x[n] = x[n] + (vx[n] + mob*fx[n]) * dt + dxb[n];
      y[n] = y[n] + (vy[n] + mob*fy[n]) * dt + dyb[n];
  
  
      // Apply boundary conditions (elastic walls)
      y[n] = std::min (h, std::max (0.0, y[n]));
      
      // Write particle position and time to the output file    
      outFile << t << "," << x[n] << "," << y[n] << "\n"; 
    
      }
      
   /*
   
   if (mod (it, Nt/100) == 1)
    figure(1)
    plot (x, y, 'ko')
    axis ([-10, 250, 0, 1])
    print ("-dpng", sprintf("frames/f_%3.3d.png", iff++))
   drawnow
  endif
endfor
*/

 /*  std::string file_temp = "frames/f_%3.3d.png"
   std::cout << file_temp << '\n'; 
   
   ofstream myfile;
   myfile.open ("example.txt");
   myfile << "Writing this to a file.\n";
   myfile.close();
   
	*/
   
   
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
