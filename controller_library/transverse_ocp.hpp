#ifndef TRANSVERSE_OCP_H_
#define TRANSVERSE_OCP_H_

#include "ocp.hpp"

extern "C" {
#include "casadi_interface.h"
}

class TransverseOcp : public Ocp {
public:
    TransverseOcp();
    virtual ~TransverseOcp();

    void Shift() override;

    double ComputeLyapunovBound(double epsN2);
    std::vector<double> ComputeTunnelBound(int Nta, double epsN2);
protected:
    casadi_memory *tunnel_bound_memory;
};

#endif // TRANSVERSE_OCP_H_