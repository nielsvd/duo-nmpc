#ifndef TANGENTIAL_OCP_H_
#define TANGENTIAL_OCP_H_

#include "ocp.hpp"

class TangentialOcp : public Ocp {
public:
    TangentialOcp();
    virtual ~TangentialOcp();

    void Shift() override;

    void SetLyapunovBound(double value);
    void SetTunnelBound(const std::vector<double> &value);

    int GetNta();
    int GetNX(int stage);
    int GetNU(int stage);
protected:
};

#endif // TANGENTIAL_OCP_H_