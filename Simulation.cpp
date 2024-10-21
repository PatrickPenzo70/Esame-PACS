#include "Simulation.h"
#include <iostream>
#include <fstream>

Simulation::Simulation(const Flow& flowProfile, double timeStep, double totalTime, double pipeRadius)
    : flow(flowProfile), timeStep(timeStep), totalTime(totalTime), pipeRadius(pipeRadius) {}

void Simulation::addBall(const Ball& ball) {
    balls.push_back(ball);
}

void Simulation::savePositionsToFile(std::ofstream& file, double time) const {
    // Write positions of all balls at the current time step
    for (size_t i = 0; i < balls.size(); ++i) {
        file << time << "," << i << ","  // Time and Ball index
             << balls[i].x << "," << balls[i].y << "," << "\n";  // Position
    }
}

void Simulation::run() {
    // Open a file to save positions
    std::ofstream file("positions.csv");
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file for writing positions.\n";
        return;
    }

    // Write header to the file
    file << "Time,BallID,PosX,PosY\n";

    double time = 0.0;
    while (time < totalTime) {
        for (auto& ball : balls) {
            ball.resetForce();  // Reset forces at the beginning of each time step

            // Apply forces
            ball.applyBrownianForce(timeStep);
            ball.applyFlowForce(flow.getVelocity(ball.y));
            ball.applyGravity();
        }

        // Handle collisions and wall bounces
        collision.handleBallCollisions(balls);
        for (auto& ball : balls) {
            collision.handleWallBounce(ball, pipeRadius);
        }

        // First half of Verlet scheme: Update positions
        for (auto& ball : balls) {
            ball.updatePosition(timeStep);
        }

        // Recompute forces based on new positions
        for (auto& ball : balls) {
            ball.applyBrownianForce(timeStep);
            ball.applyFlowForce(flow.getVelocity(ball.y));
            ball.applyGravity();
        }

        // Second half of Verlet scheme: Update velocities
        for (auto& ball : balls) {
            ball.updateVelocity(timeStep);
        }

        // Save the current positions of the balls to the file
        savePositionsToFile(file, time);

        // Advance simulation time
        time += timeStep;
    }

    // Close the file after the simulation
    file.close();
    std::cout << "Simulation complete. Positions saved to 'positions.csv'.\n";
}
