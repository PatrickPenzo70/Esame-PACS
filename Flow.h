#ifndef FLOW_H
#define FLOW_H

class Flow {
private:
    double maxVelocity;
    double pipeRadius;

public:
    Flow(double maxVelocity, double pipeRadius);

    double getVelocity(double y) const;
};

#endif // FLOW_H
