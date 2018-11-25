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

#include "transverse_ocp.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>

extern "C" {
#include "codegen/trans_dims_bounds.h"
#include "codegen/trans_functions.h"
#include <acados/ocp_qp/ocp_qp_hpipm.h>
}

TransverseOcp::TransverseOcp() : Ocp(getN(), getNP()) {
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

    dims = (sqp_dims *)malloc(sizeof(sqp_dims));
    dims->N = N;
    dims->nx = nx;
    dims->nu = nu;
    dims->nb = nb;
    dims->ng = ng;

    opts = prepare_sqp_opts(dims);
    opts->eh = true;
    opts->its = 10;
    opts->tol = 1e-6;
    opts->reg_eps = 1e-8;
    for (int i=0; i<N; ++i) {
        // jk
        opts->jk_fun[i] = &trans_jk;
        opts->jk_sparsity_in[i] = &trans_jk_sparsity_in;
        opts->jk_sparsity_out[i] = &trans_jk_sparsity_out;
        opts->jk_work[i] = &trans_jk_work;
        // gk
        opts->gk_fun[i] = &trans_gk;
        opts->gk_sparsity_in[i] = &trans_gk_sparsity_in;
        opts->gk_sparsity_out[i] = &trans_gk_sparsity_out;
        opts->gk_work[i] = &trans_gk_work;
        // fk
        opts->fk_fun[i] = &trans_fk;
        opts->fk_sparsity_in[i] = &trans_fk_sparsity_in;
        opts->fk_sparsity_out[i] = &trans_fk_sparsity_out;
        opts->fk_work[i] = &trans_fk_work;
    }
    // jN
    opts->jk_fun[N] = &trans_jN;
    opts->jk_sparsity_in[N] = &trans_jN_sparsity_in;
    opts->jk_sparsity_out[N] = &trans_jN_sparsity_out;
    opts->jk_work[N] = &trans_jN_work;
    // gN
    opts->gk_fun[N] = &trans_gN;
    opts->gk_sparsity_in[N] = &trans_gN_sparsity_in;
    opts->gk_sparsity_out[N] = &trans_gN_sparsity_out;
    opts->gk_work[N] = &trans_gN_work;

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

    tunnel_bound_memory = data->jk_memory[N];

    free(qpconfig);
}

TransverseOcp::~TransverseOcp() {
    free(qpopts);
    free(dims);
    free(opts);
    free(in);
    free(out);
    free(data);
}

void TransverseOcp::Shift() {
    for (int i=0;i<N;++i) {
        for (int ii=0;ii<data->qpdims->nx[i];++ii) out->x[i][ii] = out->x[i+1][ii];
    }
    for (int i=0;i<N-1;++i) {
        for (int ii=0;ii<data->qpdims->nu[i];++ii) out->u[i][ii] = out->u[i+1][ii];
        for (int ii=0;ii<data->qpdims->nx[i+1];++ii) out->pi[i][ii] = out->pi[i+1][ii];
    }
}

double TransverseOcp::ComputeLyapunovBound(double epsN2) {
    int N = getN();
    double ret = 0;
    for (int i=1; i<=N; ++i) ret += data->stage_cost[i];
    return fmax(ret, epsN2);
}

std::vector<double> TransverseOcp::ComputeTunnelBound(int Nta, double epsN2) {
    int N = getN();
    if ((Nta < 0) || (Nta > N)) {
        std::cout << "Nta < 0 does not make sense, Nta > Ntr is not yet supported." << std::endl;
        return std::vector<double>(0);
    }
    std::vector<double> ret(Nta, 0);

    for (int i=0; i<Nta; ++i) {
        eval_jk(out->x[i+1], out->u[i+1], data->lag_hess[i+1], in->p, NULL, NULL, NULL, ret.data()+i, data->jk_memory[N]);
        ret[i] = fmax(ret[i], epsN2);
    }

    return ret;
}