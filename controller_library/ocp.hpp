#ifndef OCP_H_
#define OCP_H_

extern "C" {
#include "sqp.h"
}

#include <vector>

class Ocp {
public:
    Ocp(int _N, int _np);
    virtual ~Ocp();

    int SolveOcp(const std::vector<double> &x0, const std::vector<double> &p, int sqp_its, double reg_eps);

    virtual void Shift();
    
    double GetKkt();

    void SetX(int stage, const std::vector<double> &value);
    void SetU(int stage, const std::vector<double> &value);
    void SetX(int stage, const std::vector<int> &idcs, const std::vector<double> &value);
    void SetU(int stage, const std::vector<int> &idcs, const std::vector<double> &value);

    std::vector<double> GetX(int stage);
    std::vector<double> GetU(int stage);
    std::vector<double> GetLamLb(int stage);
    std::vector<double> GetLamUb(int stage);
    std::vector<double> GetLamLg(int stage);
    std::vector<double> GetLamUg(int stage);
    std::vector<double> GetPi(int stage);

    int N;

  protected:
    void *qpopts;
    sqp_dims *dims;
    sqp_opts *opts;
    sqp_in *in;
    sqp_out *out;
    sqp_data *data;

    int np;
};

#endif // OCP_H_