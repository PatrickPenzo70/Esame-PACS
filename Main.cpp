#include "ConfigReader.h"
#include "Simulation.h"

int main() {
    // Load configuration
    ConfigReader config("config/config.cfg");

    // Retrieve values from the configuration
    int numBalls = config.getInt("num_balls");
    double timeStep = config.getDouble("time_step");
    double totalTime = config.getDouble("total_time");
    double pipeRadius = config.getDouble("pipe_radius");
    double maxFlowVelocity = config.getDouble("flow_max_velocity");

    // Initialize the flow and simulation
    Flow flow(maxFlowVelocity, pipeRadius);
    Simulation simulation(flow, timeStep, totalTime, pipeRadius);

    // Initialize balls (with arbitrary positions for this example)
    for (int i = 0; i < numBalls; ++i) {
        Ball ball(0.0, 0.1 * i, 0.05, 1.0);  // Example initialization
        simulation.addBall(ball);
    }

    // Run the simulation
    simulation.run();

    return 0;
}

