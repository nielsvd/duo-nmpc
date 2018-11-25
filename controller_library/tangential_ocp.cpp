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

#include "tangential_ocp.hpp"

#include <cstdlib>

extern "C" {
#include "codegen/tang_dims_bounds.h"
#include "codegen/tang_functions.h"
#include <acados/ocp_qp/ocp_qp_hpipm.h>
}

TangentialOcp::TangentialOcp() : Ocp(getN(), getNP()) {
    int N = getN();
    int tmp;
    void *ptr;
    qp_solver_config *qpconfig;
    tmp = ocp_qp_solver_config_calculate_size();
    ptr = calloc(1, tmp);
    qpconfig = ocp_qp_solver_config_assign(ptr);
    ocp_qp_hpipm_config_initialize_default(qpconfig);

    ocp_qp_dims qpdims;
    qpdims.N = N;
    qpdims.nx = qpnx;
    qpdims.nu = nu;
    qpdims.nbx = nbx;
    qpdims.nbu = nbu;
    qpdims.nb = nb;
    qpdims.ng = ng;
    qpdims.ns = ns;

    tmp = ocp_qp_hpipm_opts_calculate_size(qpconfig, &qpdims);
    ptr = calloc(1, tmp);
    qpopts = ocp_qp_hpipm_opts_assign(qpconfig, &qpdims, ptr);
    ocp_qp_hpipm_opts_initialize_default(qpconfig, &qpdims, qpopts);

    dims = (sqp_dims *) malloc(sizeof(sqp_dims));
    dims->N = N;
    dims->nx = nx;
    dims->nu = nu;
    dims->nb = nb;
    dims->ng = ng;

    opts = prepare_sqp_opts(dims);
    opts->eh = true;
    opts->its = 10;
    opts->tol = 1e-6;
    opts->reg_eps = 1e-6;
    for (int i=0; i<N; ++i) {
        // jk
        opts->jk_fun[i] = tang_jk(i);
        opts->jk_sparsity_in[i] = tang_jk_sparsity_in(i);
        opts->jk_sparsity_out[i] = tang_jk_sparsity_out(i);
        opts->jk_work[i] = tang_jk_work(i);
        // gk
        opts->gk_fun[i] = tang_gk(i);
        opts->gk_sparsity_in[i] = tang_gk_sparsity_in(i);
        opts->gk_sparsity_out[i] = tang_gk_sparsity_out(i);
        opts->gk_work[i] = tang_gk_work(i);
        // fk
        opts->fk_fun[i] = tang_fk(i);
        opts->fk_sparsity_in[i] = tang_fk_sparsity_in(i);
        opts->fk_sparsity_out[i] = tang_fk_sparsity_out(i);
        opts->fk_work[i] = tang_fk_work(i);
    }
    // jk
    opts->jk_fun[N] = tang_jk(N);
    opts->jk_sparsity_in[N] = tang_jk_sparsity_in(N);
    opts->jk_sparsity_out[N] = tang_jk_sparsity_out(N);
    opts->jk_work[N] = tang_jk_work(N);
    // gk
    opts->gk_fun[N] = tang_gk(N);
    opts->gk_sparsity_in[N] = tang_gk_sparsity_in(N);
    opts->gk_sparsity_out[N] = tang_gk_sparsity_out(N);
    opts->gk_work[N] = tang_gk_work(N);

    // sqp_in
    in = (sqp_in *) malloc(sizeof(sqp_in));
    in->idxb = idxb;
    in->lb = lb;
    in->ub = ub;
    in->lg = lg;
    in->ug = ug;

    // sqp_out
    out = prepare_sqp_out(&qpdims);
    // Zero-out initial solution
    for (int i=0;i<=N;++i) {
        for (int ii=0;ii<qpdims.nx[i];++ii) out->x[i][ii] = 0;
        for (int ii=0;ii<qpdims.nu[i];++ii) out->u[i][ii] = 0;
        for (int ii=0;ii<qpdims.nb[i];++ii) out->lam_lb[i][ii] = 0;
        for (int ii=0;ii<qpdims.nb[i];++ii) out->lam_ub[i][ii] = 0;
        for (int ii=0;ii<qpdims.ng[i];++ii) out->lam_lg[i][ii] = 0;
        for (int ii=0;ii<qpdims.ng[i];++ii) out->lam_ug[i][ii] = 0;
    }
    for (int i=0;i<N;++i) {
        for (int ii=0;ii<qpdims.nx[i+1];++ii) out->pi[i][ii] = 0;
    }

    data = prepare_sqp_data(dims, opts, &qpdims, qpconfig, qpopts);

    free(qpconfig);
}

TangentialOcp::~TangentialOcp() {
    free(qpopts);
    free(dims);
    free(opts);
    free(in);
    free(out);
    free(data);
}

void TangentialOcp::Shift() {
    int nx_basis = dims->nx[0];
    int nu_basis = dims->nu[0];
    
    for (int i=0;i<getNta();++i) {
        for (int ii=0;ii<nx_basis;++ii) out->x[i][ii] = out->x[i+1][ii];
    }
    for (int i=0;i<getNta()-1;++i) {
        for (int ii=0;ii<nu_basis;++ii) out->u[i][ii] = out->u[i+1][ii];
        for (int ii=0;ii<nx_basis;++ii) out->pi[i][ii] = out->pi[i+1][ii];
    }
}

void TangentialOcp::SetLyapunovBound(double value) {
    ugdata[lyapbound_idx] = value;
}

void TangentialOcp::SetTunnelBound(const std::vector<double> &value) {
    int Nta = getNta();
    for (int i=0; i<Nta; ++i) ugdata[tunnelbound_idx[i]] = value[i];
}

int TangentialOcp::GetNta() {
    return getNta();
}

int TangentialOcp::GetNX(int stage) {
    return dims->nx[stage];
}

int TangentialOcp::GetNU(int stage) {
    return dims->nu[stage];
}