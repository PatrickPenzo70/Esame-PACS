#include "Balls.h"
#include <cmath>    // For random Brownian force

Ball::Ball(double x_init, double y_init, double r, double m)
    : x(x_init), y(y_init), radius(r), mass(m), vx(0), vy(0), fx(0), fy(0) {}

void Ball::updatePosition(double dt) {
    // Verlet integration scheme (assuming velocity Verlet)
    x += vx * dt + 0.5 * (fx / mass) * dt * dt;
    y += vy * dt + 0.5 * (fy / mass) * dt * dt;

    // Update velocity (assuming the velocity is updated in the simulation after force update)
    vx += 0.5 * (fx / mass) * dt;
    vy += 0.5 * (fy / mass) * dt;
}

void Ball::applyBrownianForce(double dt) {
    // Brownian force (simple model using random force)
    double noise_strength = 0.1;  // Tunable parameter
    fx += noise_strength * ((rand() / (double)RAND_MAX) - 0.5);
    fy += noise_strength * ((rand() / (double)RAND_MAX) - 0.5);
}

void Ball::applyGravity() {
    // Apply the gravitational force in the negative z-direction
    const double gravity = 9.81;  // Gravitational acceleration (m/s^2)
    fy -= mass * gravity;         // Force due to gravity is mass * g in the negative y-direction

}

void Ball::applyFlowForce(double flowVelocity) {
    // Drag force approximation based on Poiseuille flow
    double dragCoefficient = 6 * M_PI * 0.01 * radius;  // Example with a constant drag coefficient
    fx += -dragCoefficient * (vx - flowVelocity);  // Fluid velocity affects the force
}

void Ball::resetForce() {
    fx = 0;
    fy = 0;
}
