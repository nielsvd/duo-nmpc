#include "ocp.hpp"

#include <cstdlib>
#include <cstring>
#include <iostream>

Ocp::Ocp(int _N, int _np) : N(_N), np(_np)  {

}

Ocp::~Ocp() {

}

int Ocp::SolveOcp(const std::vector<double> &x0, const std::vector<double> &p, int sqp_its, double reg_eps) {
    // Set sqp_in
    in->x0 = x0.data();

    if (p.size() < np) {
        std::cout << "SolveOcp: dimension mismatch." << std::endl;
        return -1;
    }

    in->p = p.data();

    // Set no. SQP iterations
    opts->its = sqp_its;
    opts->reg_eps = reg_eps;

    // Solve OCP
    int status = sqp(dims, in, out, opts, data);

    return status;
}

void Ocp::Shift() {
    
}

double Ocp::GetKkt() {
    return compute_kkt(data);
}

void Ocp::SetX(int stage, const std::vector<double>& value) {
    if ((stage < 0) || (stage > N)) {
        std::cout << "SetX: stage out of bounds." << std::endl;
        return;
    }
    int size = data->qpdims->nx[stage];
    if (value.size() != size) {
        std::cout << "SetX: dimension mismatch." << std::endl;
        return;      
    }
    memcpy(out->x[stage],value.data(),size*sizeof(double));
}

void Ocp::SetU(int stage, const std::vector<double>& value) {
    if ((stage < 0) || (stage >= N)) {
        std::cout << "SetU: stage out of bounds." << std::endl;
        return;
    }
    int size = data->qpdims->nu[stage];
    if (value.size() != size) {
        std::cout << "SetU: dimension mismatch." << std::endl;
        return;      
    }
    memcpy(out->u[stage],value.data(),size*sizeof(double));
}

void Ocp::SetX(int stage, const std::vector<int> &idcs, const std::vector<double> &value) {
    if ((stage < 0) || (stage > N)) {
        std::cout << "SetX: stage out of bounds." << std::endl;
        return;
    }
    int size = idcs.size();
    for (auto idx : idcs) {
        if (idx >= data->qpdims->nx[stage]) {
            std::cout << "SetX: dimension mismatch." << std::endl;
            return;      
        }
    }
    for (int i=0; i<size; ++i) out->x[stage][idcs[i]] = value[i];
}

void Ocp::SetU(int stage, const std::vector<int> &idcs, const std::vector<double> &value) {
    if ((stage < 0) || (stage >= N)) {
        std::cout << "SetU: stage out of bounds." << std::endl;
        return;
    }
    int size = idcs.size();
    for (auto idx : idcs) {
        if (idx >= data->qpdims->nu[stage]) {
            std::cout << "SetU: dimension mismatch." << std::endl;
            return;      
        }
    }
    for (int i=0; i<size; ++i) out->u[stage][idcs[i]] = value[i];
}

std::vector<double> Ocp::GetX(int stage) {
    if ((stage < 0) || (stage > N)) {
        std::cout << "GetX: stage out of bounds." << std::endl;
        return std::vector<double>(0);
    }

    int size = data->qpdims->nx[stage];
    std::vector<double> ret(size);
    memcpy(ret.data(),out->x[stage],size*sizeof(double));
    return ret;
}

std::vector<double> Ocp::GetU(int stage) {
    if ((stage < 0) || (stage > N)) {
        std::cout << "GetU: stage out of bounds." << std::endl;
        return std::vector<double>(0);
    }

    int size = data->qpdims->nu[stage];
    std::vector<double> ret(size);
    memcpy(ret.data(),out->u[stage],size*sizeof(double));
    return ret;
}

std::vector<double> Ocp::GetLamLb(int stage) {
    if ((stage < 0) || (stage > N)) {
        std::cout << "GetLamLb: stage out of bounds." << std::endl;
        return std::vector<double>(0);
    }

    int size = data->qpdims->nb[stage];
    std::vector<double> ret(size);
    memcpy(ret.data(),out->lam_lb[stage],size*sizeof(double));
    return ret;
}

std::vector<double> Ocp::GetLamUb(int stage) {
    if ((stage < 0) || (stage > N)) {
        std::cout << "GetLamUb: stage out of bounds." << std::endl;
        return std::vector<double>(0);
    }

    int size = data->qpdims->nb[stage];
    std::vector<double> ret(size);
    memcpy(ret.data(),out->lam_ub[stage],size*sizeof(double));
    return ret;
}

std::vector<double> Ocp::GetLamLg(int stage) {
    if ((stage < 0) || (stage > N)) {
        std::cout << "GetLamLg: stage out of bounds." << std::endl;
        return std::vector<double>(0);
    }

    int size = data->qpdims->ng[stage];
    std::vector<double> ret(size);
    memcpy(ret.data(),out->lam_lg[stage],size*sizeof(double));
    return ret;
}

std::vector<double> Ocp::GetLamUg(int stage) {
    if ((stage < 0) || (stage > N)) {
        std::cout << "GetLamUg: stage out of bounds." << std::endl;
        return std::vector<double>(0);
    }

    int size = data->qpdims->ng[stage];
    std::vector<double> ret(size);
    memcpy(ret.data(),out->lam_ug[stage],size*sizeof(double));
    return ret;
}

std::vector<double> Ocp::GetPi(int stage) {
    if ((stage < 0) || (stage > N-1)) {
        std::cout << "GetLamUg: stage out of bounds." << std::endl;
        return std::vector<double>(0);
    }

    int size = data->qpdims->nx[stage + 1];
    std::vector<double> ret(size);
    memcpy(ret.data(),out->pi[stage],size*sizeof(double));
    return ret;
}