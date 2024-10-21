#include "Flow.h"

Flow::Flow(double maxVelocity, double pipeRadius)
    : maxVelocity(maxVelocity), pipeRadius(pipeRadius) {}

double Flow::getVelocity(double y) const {
    // Parabolic velocity profile (Poiseuille flow)
    return maxVelocity * (1 - (y * y) / (pipeRadius * pipeRadius));
}
