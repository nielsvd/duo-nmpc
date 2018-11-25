/*
 *    This file is part of duo-nmpc.
 *
 *    Duo-nmpc
 *    Copyright (C) 2018 Niels van Duijkeren, KU Leuven.
 *
 *    Duo-nmpc is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    Duo-nmpc is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with duo-nmpc; if not, write to the Free Software Foundation,
 *    Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#include "duo_nmpc.hpp"

DuoNmpc::DuoNmpc(const std::vector<double> &xinit, const std::vector<double> &uinit, int maxit, double reg_eps) : maxit(maxit), reg_eps(reg_eps), reset_tang(false) {
    sys_nx = xinit.size();
    sys_nu = uinit.size();

    for (int i=1; i<=ocp_tr.N; ++i) ocp_tr.SetX(i,xinit);
    for (int i=0; i<ocp_tr.N; ++i) ocp_tr.SetU(i,uinit);

    // Assuming Nta <= Ntr for now!!!!
    for (int i=1; i<=ocp_ta.N; ++i) {
        std::vector<double> x0i(ocp_ta.GetNX(i),0);
        if (i <= ocp_ta.GetNta()) {
            for (int ii=0; ii<xinit.size(); ++ii) x0i[ii] = xinit[ii];
        }
        ocp_ta.SetX(i, x0i);
    };

    for (int i=0; i<ocp_ta.N; ++i) {
        std::vector<double> u0i(ocp_ta.GetNU(i),0);
        if (i < ocp_ta.GetNta() && i > 0) {
            for (int ii=0; ii<uinit.size(); ++ii) u0i[ii] = uinit[ii];
        }
        ocp_ta.SetU(i, u0i);
    };
}

DuoNmpc::~DuoNmpc() {

}

std::vector<double> DuoNmpc::Solve(const std::vector<double> &x0, const std::vector<double> &ptr, const std::vector<double> &pta, double epsN) {
    // Solve transverse problem
    int status_tr = ocp_tr.SolveOcp(x0, ptr, maxit, reg_eps);

    // Compute Tunnel and Lyapunov bounds
    double epsN2 = epsN*epsN;
    std::vector<double> tunnel_bound = ocp_tr.ComputeTunnelBound(ocp_ta.GetNta(), epsN2);
    double lyap_bound = ocp_tr.ComputeLyapunovBound(epsN2);

    // Set Tunnel and Lyapunov bounds
    ocp_ta.SetLyapunovBound(lyap_bound);
    ocp_ta.SetTunnelBound(tunnel_bound);

    // Set solution of transverse OCP as initial guess to tangential OCP
    //      Assuming Nta <= Ntr for now!!!!
    for (int i=1; i<=ocp_tr.N; ++i) {
        std::vector<int> idcs(sys_nx,0);
        int offset = (i <= ocp_ta.GetNta()) ? sys_nx : 0;
        for (int ii=0; ii<sys_nx; ++ii) idcs[ii] = offset + ii;
        ocp_ta.SetX(i, idcs, ocp_tr.GetX(i));
    };

    for (int i=1; i<ocp_tr.N; ++i) {
        std::vector<int> idcs(sys_nu,0);
        int offset = (i < ocp_ta.GetNta()) ? sys_nu : 0;
        for (int ii=0; ii<sys_nu; ++ii) idcs[ii] = offset + ii;
        ocp_ta.SetU(i, idcs, ocp_tr.GetU(i));
    };

    if (reset_tang) {
        reset_tang = false;
        for (int i=1; i<=ocp_ta.N; ++i) {
            std::vector<int> idcs(sys_nx,0);
            for (int ii=0; ii<sys_nx; ++ii) idcs[ii] = ii;
            ocp_ta.SetX(i, idcs, ocp_tr.GetX(i));
        };
        for (int i=0; i<ocp_ta.N; ++i) {
            std::vector<int> idcs(sys_nu,0);
            for (int ii=0; ii<sys_nu; ++ii) idcs[ii] = ii;
            ocp_ta.SetU(i, idcs, ocp_tr.GetU(i));
        };
    }

    // Solve tangential problem
    int status_ta = ocp_ta.SolveOcp(x0, pta, maxit, reg_eps);

    std::vector<double> u;
    // Select control
    if (status_ta==0) {
        u = ocp_ta.GetU(0);
    } else {
        u = ocp_tr.GetU(0);
        reset_tang = true;
    }

    // return economic control
    return u;
}

void DuoNmpc::Shift() {
    ocp_tr.Shift();
    ocp_ta.Shift();
}

double DuoNmpc::GetTransverseKkt() {
    return ocp_tr.GetKkt();
}

double DuoNmpc::GetTangentialKkt() {
    return ocp_ta.GetKkt();
}
