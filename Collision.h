#ifndef COLLISION_H
#define COLLISION_H

#include "Balls.h"
#include <vector>

class Collision {
public:
    void handleBallCollisions(std::vector<Ball>& balls);
    void handleWallBounce(Ball& ball, double pipeRadius);
};

#endif // COLLISION_H
