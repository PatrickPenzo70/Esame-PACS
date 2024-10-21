#include "Collision.h"
#include <cmath>

void Collision::handleBallCollisions(std::vector<Ball>& balls) {
    for (size_t i = 0; i < balls.size(); ++i) {
        for (size_t j = i + 1; j < balls.size(); ++j) {
            double dx = balls[j].x - balls[i].x;
            double dy = balls[j].y - balls[i].y;
            double dist = std::sqrt(dx * dx + dy * dy);
            double minDist = balls[i].radius + balls[j].radius;

            if (dist < minDist) {
                // Elastic collision response between balls
                double overlap = minDist - dist;
                double nx = dx / dist, ny = dy / dist;

                // Simple response: push balls apart and invert velocity
                balls[i].vx -= overlap * nx;
                balls[i].vy -= overlap * ny;
                balls[j].vx += overlap * nx;
                balls[j].vy += overlap * ny;
            }
        }
    }
}

void Collision::handleWallBounce(Ball& ball, double pipeRadius) {
    // Check against the pipe walls (bounce back if beyond pipe radius)
    if (std::fabs(ball.y) > pipeRadius - ball.radius) {
        ball.vy = -ball.vy;  // Invert velocity in y-direction
    }
}
