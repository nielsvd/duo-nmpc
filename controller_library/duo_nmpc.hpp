#ifndef DUO_NMPC_HPP_
#define DUO_NMPC_HPP_

#include <vector>
#include "transverse_ocp.hpp"
#include "tangential_ocp.hpp"

class DuoNmpc {
public:
    DuoNmpc(const std::vector<double> &xinit, const std::vector<double> &uinit, int maxit, double reg_eps);
    virtual ~DuoNmpc();

    std::vector<double> Solve(const std::vector<double> &x0, const std::vector<double> &ptr, const std::vector<double> &pta, double epsN);
    void Shift();

    double GetTransverseKkt();
    double GetTangentialKkt();

    TransverseOcp ocp_tr;
    TangentialOcp ocp_ta;
private:
    int sys_nx;
    int sys_nu;
    int maxit;
    double reg_eps;

    bool reset_tang;
};

#endif // DUO_NMPC_HPP_