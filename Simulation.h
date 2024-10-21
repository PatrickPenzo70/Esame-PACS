#ifndef SIMULATION_H
#define SIMULATION_H

#include <vector>
#include "Balls.h"
#include "Flow.h"
#include "Collision.h"

class Simulation {
public:
    Simulation(const Flow& flowProfile, double timeStep, double totalTime, double pipeRadius);

    void addBall(const Ball& ball);
    void run();
    void savePositionsToFile(std::ofstream& file, double time) const; 

private:
    std::vector<Ball> balls;
    Flow flow;
    Collision collision;
    double timeStep;
    double totalTime;
    double pipeRadius;
};

#endif // SIMULATION_H

