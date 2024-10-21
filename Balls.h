#ifndef BALL_H
#define BALL_H

class Ball {
public:
    double x, y;   // Position of the ball
    double vx, vy; // Velocity
    double fx, fy; // Force
    double mass, radius;

    Ball(double x_init, double y_init, double r, double m);

    void applyGravity();
    void updatePosition(double dt);  // Verlet scheme
    void applyBrownianForce(double dt);
    void applyFlowForce(double flowVelocity);  // Poiseuille Flow effect
    void resetForce();  // Clear force for next step
    void updateVelocity(double timestep);
};

#endif // BALL_H

