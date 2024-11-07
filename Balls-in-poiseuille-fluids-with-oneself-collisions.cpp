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


//Constants
const double U_max = 1.0;       // Maximum fluid velcosity
const double R = 1.0;           // Pipe radius
const double mu = 0.01;         // Dinamic viscosity of the fluid
const double g = 9.81;          // Gravitational acceleration
const double dt = 0.01;         // Time step
const double totalTime = 10.0;  // Totale simulation time
const double N = 100;           // Number of particles
const double k_B = 1.38e-23;	// Boltzmann constant (J/K)
const double T = 300.0;		// Temperature (K)


// Random number generation
std::random_device rd;
std::mt19937 generator(rd()); // Mersenne Twister RNG
std::normal_distribution<double> distribution(0.0, 1.0): // Normal distribution with mean 0 and variance 1

// Ball properties
struct Ball {
    double x, y;            // Position of the particle
    double vx, vy;          // Velocity of the particle
    double r;               // Radius of the particle
    double m;               // Mass of the particle
};

// Poiseuille flow velocity at a given y position
double fluidVelocity(double y) {
    return U_max * (1 - (y * y) / R*R);
}

// Drag force at a given y position on the ball
double dragForce(double v_fluid, double v_ball, double r) {
    return 6 * M_PI * mu * r * (v_fluid - v_ball);
}

// Distance between two particles
double distance(const Ball &ball1, const Ball &ball2) {
    return sqrt((ball1.x - ball2.x) * (ball1.x - ball2.x) +
                (ball1.y - ball2.y) * (ball1.y - ball2.y));
}

// Check if two particles are colliding
bool checkCollision(const Ball &ball1, const Ball &ball2) {
    return distance(ball1, ball2) <= (ball1.r + ball2.r);
}

// Handle particle-particle elastic collision
void handleCollision(Ball &ball1, Ball &ball2) {

    // Vector from particle1 to particle2
    double dx = ball2.x - ball1.x;
    double dy = ball2.y - ball1.y;
    
    double dist = sqrt(dx*dx + dy*dy);

    // Normal vector
    double nx = dx / dist;
    double ny = dy / dist;
    
    // Relative velocity
    double dvx = ball1.vx - ball2.vx;
    double dvy = ball1.vy - ball2.vy;
    
    // Dot product of relative velocity and normal vector
    double dot = dvx * nx + dvy * ny;
    
    // Only proceed if balls are moving towards each other
    if (dot > 0) return;
    
    // Compute impulse scalar
    double impulse = (2 * dot) / (ball1.m + ball2.m);
    
    // Update velocities
    ball1.vx -= impulse * ball2.m * nx;
    ball1.vy -= impulse * ball2.m * ny;
    
    ball2.vx += impulse * ball1.m * nx;
    ball2.vy += impulse * ball1.m * ny;
    
    // Adjust positions to avoid overlap
    double overlap = 0.5 * (ball1.r + ball2.r - dist);
    ball1.x -= overlap * nx;
    ball1.y -= overlap * ny;
    
    ball2.x += overlap * nx;
    ball2.y += overlap * ny;
}

// Generate a Brownian force for a ball
void addBrownianForce(Ball &ball) {
	// Diffusion coefficient D using Stokes-Einstein relation
	double D =(k_B * T) / (6 * M_PI * mu * ball.r); 

	// Standard deviation for Brownian force
	double sigma = sqrt(2 * D / dt);

	// Generate random Brownian forces (Gaussian distributed)
	double Fx = sigma * distribution(generator);
	double Fy = sigma * distribution(generator);

	// Update velocities with Brownian force
	ball.vx += Fx / ball.m;
	ball.vy += Fy / ball.m;
}


// Update the position and velocity of the particle system
void updateBall(Ball &ball, double dt) {
    // Get fluid velocity at the ball's position
    double v_fluid = fluidVelocity(ball.y);

    // Calculate drag forces in each direction (assuming drag only in the x direction)
    double F_drag_x = dragForce(v_fluid, ball.vx, ball.r); 
    double F_drag_y = -dragForce(0, ball.vy, ball.r); // Drag due to vertical motion

    // Newton's second law to update the velocity
    ball.vx -= (F_drag_x / ball.m) * dt;
    ball.vy -= (F_drag_y / ball.m - g) * dt; // Include gravity

    // Update the position
    ball.x += ball.vx * dt;
    ball.y += ball.vy * dt;    
    
    // Bounce back if the ball hits the top or bottom of the pipe
    if (ball.y >= R) {
        // Ball hits the top wall, reverse y-velocity and move it inside the pipe with inelastic collisions
        ball.vy = -0.9*ball.vy;
        ball.y = -R - (ball.y - R); // Adjust position to avoid getting stuck
    }  else if (ball.y <= -R) {
        // Ball hits the bottom wall, reverse y-velocity and move it inside the pipe with inelastic collisions
        ball.vy = -0.9*ball.vy;
        ball.y = - R - (ball.y + R); // Adjust position to avoid getting stuck
    }
}

int main() {
    // Initialize the balls
    std::vector<Ball> balls(N);

    // Randomly initialize balls with different initial postions and velcoities
    for (int i = 0; i < N; ++i) {
        balls[i].x = 0.0;
        balls[i].y = 0.5 -i * 0.1; // Space the balls out along the y-axis
        balls[i].vx = 0.0;
        balls[i].vy = 0.0;
        balls[i].r = 0.1; // Radius of the ball
        balls[i].m = 1.0; // Mass of the ball
    }

    // Open a file to save the positions
    std::ofstream outFile("balls_position_collision_brownian.csv");

    // Check if the file is open
    if (!outFile.is_open()) {
        std::cerr << "Error opening output file." << std::endl;
        return 1;
    }

    // Write the header for the CSV file
    outFile << "Time";
    for (int i = 0; i < N; ++i) {
        outFile << ", Ball_" << i + 1 << "_X, Ball" << i + 1 << "_Y, Ball";
    }
    outFile << std::endl;

    // Simulate over time
    double time = 0.0;
    while (time < totalTime) {
        // Update the position for all balls
        for (int i = 0; i < N; ++i ) {
            updateBall(balls[i], dt);
        }

	// Check for collision between balls
	for (int i = 0; i < N; ++i) {
	    for (int j = i + 1; j = N; ++j) {
                if (checkCollision(balls[i], balls[j])) {
		    handleCollision(balls[i], balls[j]);
		}
	    }
	}


        // Write the current time and positions of all balls to the CSV file
        outFile << time;
        for (int i = 0; i < N; ++i) {
            outFile << ", " << balls[i].x << ", " << balls[i].y; 
        }
        outFile << std::endl;

        // Increase the time
        time += dt;
    }

    // Close the file
    outFile.close();
    
    std::cout << "Simulation completed. Particles positions saved to balls_position_collision.csv" << std::endl;
    
    return 0;
}



// DEFINITIONS:

// fluidVelocity(): Defines the velocity of the fluid as a function of the particle's y position, using
// the parabolic profile for Poiseuille flow.

// dragForce(): Calculates the drag force based on the difference  between the fluid's velocity and the
// particle's velocity.

// updateBall(): Updates the particle's position and velocity using the drag force and gravity, and
// advances the position by one time step.

// Collision Detection: The collision is detected by checking if the y-position of the ball exceeds 
// the pipe's boundaries, i.e., if y ≥ R (top wall) or y ≤ −R (bottom wall).

// Bounce-back Mechanism: When the balls hit a wall, I reverse its velocity in the y-direction: 
// ball.vy = -ball.vy. After reversing the velocity, I adjust the position so that the ball is placed 
// slightly inside the wall. This adjustment prevents the balls from being stuck inside the wall 
// due to numerical errors: 

// ball.y = R - (ball.y - R);  // For top wall
// ball.y = -R - (ball.y + R); // For bottom wall

