#include "sqp.h"

#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <acados/ocp_qp/ocp_qp_common.h>

// Debug
#define PRINTCOST 0
#define PRINTDYN 0
#define PRINTCON 0
#define PRINTSOL 0
#define CAPCONTROLS 1

sqp_dims *prepare_sqp_dims(int N) {
    //
    // Compute size of workspace
    //
    int size = sizeof(sqp_dims);
    // nx
    size += (N+1)*sizeof(int);
    // nu
    size += (N+1)*sizeof(int);
    // ng
    size += (N+1)*sizeof(int);
    // nb
    size += (N+1)*sizeof(int);

    //
    // Create workspace
    //
    char *raw_mem = (char *) malloc(size);

    //
    // Assign workspace
    //
    char *c_ptr = raw_mem;

    sqp_dims *out = (sqp_dims *)c_ptr;
    c_ptr += sizeof(sqp_dims);
    // nx
    out->nx = (int *) c_ptr;
    c_ptr += (N+1)*sizeof(int);
    // nu
    out->nu = (int *) c_ptr;
    c_ptr += (N+1)*sizeof(int);
    // ng
    out->ng = (int *) c_ptr;
    c_ptr += (N+1)*sizeof(int);
    // nb
    out->nb = (int *) c_ptr;
    c_ptr += (N+1)*sizeof(int);


    assert(c_ptr == raw_mem + size && "Assignment error.");

    // Configure
    out->N = N;

    return out;
}

sqp_opts *prepare_sqp_opts(sqp_dims *dims) {
    int N = dims->N;

    //
    // Compute size of workspace
    //
    int size = sizeof(sqp_opts);
    int tmp_size = sizeof(int (*)(int *, int *, int *, int *)); // work
    tmp_size += sizeof(const int *(*)(int)); // sparsity_in
    tmp_size += sizeof(const int *(*)(int)); // sparsity_out
    tmp_size += sizeof(int (*)(const double **, double **, int *, double *, void *)); // fun
    size += 2*(N+1)*tmp_size + N*tmp_size;

    //
    // Create workspace
    //
    char *raw_mem = (char *) malloc(size);

    //
    // Assign workspace
    //
    char *c_ptr = raw_mem;

    sqp_opts *out = (sqp_opts *)c_ptr;
    c_ptr += sizeof(sqp_opts);
    // jk
    out->jk_work = (int (**)(int *, int *, int *, int *)) c_ptr;
    c_ptr += (N+1)*sizeof(int (*)(int *, int *, int *, int *));
    out->jk_sparsity_in = (const int *(**)(int)) c_ptr;
    c_ptr += (N+1)*sizeof(const int *(*)(int));
    out->jk_sparsity_out = (const int *(**)(int)) c_ptr;
    c_ptr += (N+1)*sizeof(const int *(*)(int));
    out->jk_fun = (int (**)(const double **, double **, int *, double *, void *)) c_ptr;
    c_ptr += (N+1)*sizeof(int (*)(const double **, double **, int *, double *, void *));
    // gk
    out->gk_work = (int (**)(int *, int *, int *, int *)) c_ptr;
    c_ptr += (N+1)*sizeof(int (*)(int *, int *, int *, int *));
    out->gk_sparsity_in = (const int *(**)(int)) c_ptr;
    c_ptr += (N+1)*sizeof(const int *(*)(int));
    out->gk_sparsity_out = (const int *(**)(int)) c_ptr;
    c_ptr += (N+1)*sizeof(const int *(*)(int));
    out->gk_fun = (int (**)(const double **, double **, int *, double *, void *)) c_ptr;
    c_ptr += (N+1)*sizeof(int (*)(const double **, double **, int *, double *, void *));
    // fk
    out->fk_work = (int (**)(int *, int *, int *, int *)) c_ptr;
    c_ptr += N*sizeof(int (*)(int *, int *, int *, int *));
    out->fk_sparsity_in = (const int *(**)(int)) c_ptr;
    c_ptr += N*sizeof(const int *(*)(int));
    out->fk_sparsity_out = (const int *(**)(int)) c_ptr;
    c_ptr += N*sizeof(const int *(*)(int));
    out->fk_fun = (int (**)(const double **, double **, int *, double *, void *)) c_ptr;
    c_ptr += N*sizeof(int (*)(const double **, double **, int *, double *, void *));

    //
    // Check
    //
    assert(c_ptr == raw_mem + size && "Assignment error.");

    return out;
}

int qp_sqp_out_calculate_size(ocp_qp_dims *qpdims) {
    int N = qpdims->N;
    //
    // Compute size of workspace
    //
    int size = sizeof(qp_colmaj_out);
    size += 6*(N+1)*sizeof(double *);
    size += N*sizeof(double *);
    for (int i=0;i<=N;++i) {
        // x
        size += qpdims->nx[i]*sizeof(double);
        // u
        size += qpdims->nu[i]*sizeof(double);
        // lam_lb
        size += qpdims->nb[i]*sizeof(double);
        // lam_ub
        size += qpdims->nb[i]*sizeof(double);
        // lam_lg
        size += qpdims->ng[i]*sizeof(double);
        // lam_ug
        size += qpdims->ng[i]*sizeof(double);
    }
    for (int i=0;i<N;++i) {
        // pi
        size += qpdims->nx[i+1]*sizeof(double);
    }

    return size;
}

void *qp_sqp_out_assign(ocp_qp_dims *qpdims, void *raw_memory) {
    int N = qpdims->N;
    //
    // Assign workspace
    //
    char *c_ptr = (char *) raw_memory;

    qp_colmaj_out *out = (qp_colmaj_out *) c_ptr;
    c_ptr += sizeof(qp_colmaj_out);
    // x
    out->x = (double **) c_ptr;
    c_ptr += (N+1)*sizeof(double *);
    for (int i=0;i<=N;++i) {
        out->x[i] = (double *) c_ptr;
        c_ptr += qpdims->nx[i]*sizeof(double);
    }
    // u
    out->u = (double **) c_ptr;
    c_ptr += (N+1)*sizeof(double *);
    for (int i=0;i<=N;++i) {
        out->u[i] = (double *) c_ptr;
        c_ptr += qpdims->nu[i]*sizeof(double);
    }
    // lam_lb
    out->lam_lb = (double **) c_ptr;
    c_ptr += (N+1)*sizeof(double *);
    for (int i=0;i<=N;++i) {
        out->lam_lb[i] = (double *) c_ptr;
        c_ptr += qpdims->nb[i]*sizeof(double);
    }
    // lam_ub
    out->lam_ub = (double **) c_ptr;
    c_ptr += (N+1)*sizeof(double *);
    for (int i=0;i<=N;++i) {
        out->lam_ub[i] = (double *) c_ptr;
        c_ptr += qpdims->nb[i]*sizeof(double);
    }
    // lam_lg
    out->lam_lg = (double **) c_ptr;
    c_ptr += (N+1)*sizeof(double *);
    for (int i=0;i<=N;++i) {
        out->lam_lg[i] = (double *) c_ptr;
        c_ptr += qpdims->ng[i]*sizeof(double);
    }
    // lam_ug
    out->lam_ug = (double **) c_ptr;
    c_ptr += (N+1)*sizeof(double *);
    for (int i=0;i<=N;++i) {
        out->lam_ug[i] = (double *) c_ptr;
        c_ptr += qpdims->ng[i]*sizeof(double);
    }
    // pi
    out->pi = (double **) c_ptr;
    c_ptr += N*sizeof(double *);
    for (int i=0;i<N;++i) {
        out->pi[i] = (double *) c_ptr;
        c_ptr += qpdims->nx[i+1]*sizeof(double);
    }

    //
    // Check
    //
    assert(c_ptr == raw_memory + qp_sqp_out_calculate_size(qpdims) && "Assignment error.");

    return out;
}

void *prepare_qp_sqp_out(ocp_qp_dims *qpdims) {
    //
    // Compute size of workspace
    //
    int size = qp_sqp_out_calculate_size(qpdims);

    //
    // Create workspace
    //
    void *raw_mem = malloc(size);

    //
    // Assign workspace
    //
    void *out = qp_sqp_out_assign(qpdims, raw_mem);

    return out;
}

int qp_colmaj_in_calculate_size(ocp_qp_dims *qpdims) {
    int N = qpdims->N;
    //
    // Compute size of workspace
    //
    int size = sizeof(qp_colmaj_in);
    size += 3*N*sizeof(double *);
    size += 11*(N+1)*sizeof(double *);
    size += (N+1)*sizeof(int *);
    for (int i=0;i<N;++i) {
        // A
        size += qpdims->nx[i+1]*qpdims->nx[i]*sizeof(double);
        // B
        size += qpdims->nx[i+1]*qpdims->nu[i]*sizeof(double);
        // b
        size += qpdims->nx[i+1]*sizeof(double);
    }
    for (int i=0;i<=N;++i) {
        // Q
        size += qpdims->nx[i]*qpdims->nx[i]*sizeof(double);
        // S
        size += qpdims->nu[i]*qpdims->nx[i]*sizeof(double);
        // R
        size += qpdims->nu[i]*qpdims->nu[i]*sizeof(double);
        // q
        size += qpdims->nx[i]*sizeof(double);
        // r
        size += qpdims->nu[i]*sizeof(double);
        // idxb
        size += qpdims->nb[i]*sizeof(int);
        // lb
        size += qpdims->nb[i]*sizeof(double);
        // ub
        size += qpdims->nb[i]*sizeof(double);
        // C
        size += qpdims->ng[i]*qpdims->nx[i]*sizeof(double);
        // D
        size += qpdims->ng[i]*qpdims->nu[i]*sizeof(double);
        // lg
        size += qpdims->ng[i]*sizeof(double);
        // ug
        size += qpdims->ng[i]*sizeof(double);
    }

    return size;
}

qp_colmaj_in *qp_colmaj_in_assign(ocp_qp_dims *qpdims, void *raw_memory) {
    int N = qpdims->N;
    //
    // Assign workspace
    //
    char *c_ptr = (char *) raw_memory;

    qp_colmaj_in *out = (qp_colmaj_in *) c_ptr;
    c_ptr += sizeof(qp_colmaj_in);

    // A
    out->A = (double **) c_ptr;
    c_ptr += N*sizeof(double *);
    for (int i=0;i<N;++i) {
        out->A[i] = (double *) c_ptr;
        c_ptr += qpdims->nx[i+1]*qpdims->nx[i]*sizeof(double);
    }
    // B
    out->B = (double **) c_ptr;
    c_ptr += N*sizeof(double *);
    for (int i=0;i<N;++i) {
        out->B[i] = (double *) c_ptr;
        c_ptr += qpdims->nx[i+1]*qpdims->nu[i]*sizeof(double);
    }
    // b
    out->b = (double **) c_ptr;
    c_ptr += N*sizeof(double *);
    for (int i=0;i<N;++i) {
        out->b[i] = (double *) c_ptr;
        c_ptr += qpdims->nx[i+1]*sizeof(double);
    }
    // Q
    out->Q = (double **) c_ptr;
    c_ptr += (N+1)*sizeof(double *);
    for (int i=0;i<=N;++i) {
        out->Q[i] = (double *) c_ptr;
        c_ptr += qpdims->nx[i]*qpdims->nx[i]*sizeof(double);
    }
    // S
    out->S = (double **) c_ptr;
    c_ptr += (N+1)*sizeof(double *);
    for (int i=0;i<=N;++i) {
        out->S[i] = (double *) c_ptr;
        c_ptr += qpdims->nu[i]*qpdims->nx[i]*sizeof(double);
    }
    // R
    out->R = (double **) c_ptr;
    c_ptr += (N+1)*sizeof(double *);
    for (int i=0;i<=N;++i) {
        out->R[i] = (double *)c_ptr;
        c_ptr += qpdims->nu[i]*qpdims->nu[i]*sizeof(double);
    }
    // q
    out->q = (double **) c_ptr;
    c_ptr += (N+1)*sizeof(double *);
    for (int i=0;i<=N;++i) {
        out->q[i] = (double *) c_ptr;
        c_ptr += qpdims->nx[i]*sizeof(double);
    }
    // r
    out->r = (double **) c_ptr;
    c_ptr += (N+1)*sizeof(double *);
    for (int i=0;i<=N;++i) {
        out->r[i] = (double *) c_ptr;
        c_ptr += qpdims->nu[i]*sizeof(double);
    }
    // lb
    out->lb = (double **) c_ptr;
    c_ptr += (N+1)*sizeof(double *);
    for (int i=0;i<=N;++i) {
        out->lb[i] = (double *) c_ptr;
        c_ptr += qpdims->nb[i]*sizeof(double);
    }
    // ub
    out->ub = (double **) c_ptr;
    c_ptr += (N+1)*sizeof(double *);
    for (int i=0;i<=N;++i) {
        out->ub[i] = (double *) c_ptr;
        c_ptr += qpdims->nb[i]*sizeof(double);
    }
    // C
    out->C = (double **) c_ptr;
    c_ptr += (N+1)*sizeof(double *);
    for (int i=0;i<=N;++i) {
        out->C[i] = (double *) c_ptr;
        c_ptr += qpdims->ng[i]*qpdims->nx[i]*sizeof(double);
    }
    // D
    out->D = (double **) c_ptr;
    c_ptr += (N+1)*sizeof(double *);
    for (int i=0;i<=N;++i) {
        out->D[i] = (double *) c_ptr;
        c_ptr += qpdims->ng[i]*qpdims->nu[i]*sizeof(double);
    }
    // lg
    out->lg = (double **) c_ptr;
    c_ptr += (N+1)*sizeof(double *);
    for (int i=0;i<=N;++i) {
        out->lg[i] = (double *) c_ptr;
        c_ptr += qpdims->ng[i]*sizeof(double);
    }
    // ug
    out->ug = (double **) c_ptr;
    c_ptr += (N+1)*sizeof(double *);
    for (int i=0;i<=N;++i) {
        out->ug[i] = (double *) c_ptr;
        c_ptr += qpdims->ng[i]*sizeof(double);
    }
    // idxb (assign last to avoid alignment issues)
    out->idxb = (int **) c_ptr;
    c_ptr += (N+1)*sizeof(int *);
    for (int i=0;i<=N;++i) {
        out->idxb[i] = (int *) c_ptr;
        c_ptr += qpdims->nb[i]*sizeof(int);
    }

    //
    // Check
    //
    assert(c_ptr == raw_memory + qp_colmaj_in_calculate_size(qpdims) && "Assignment error.");

    return out;
}

qp_colmaj_in *prepare_qp_colmaj_in(ocp_qp_dims *qpdims) {
    int N = qpdims->N;
    //
    // Compute size of workspace
    //
    int size = qp_colmaj_in_calculate_size(qpdims);

    //
    // Create workspace
    //
    char *raw_mem = (char *) malloc(size);

    //
    // Assign workspace
    //
    qp_colmaj_in *out = qp_colmaj_in_assign(qpdims, raw_mem);

    return out;
}

int qp_colmaj_out_calculate_size(ocp_qp_dims *qpdims) {
    return qp_sqp_out_calculate_size(qpdims);
}

qp_colmaj_out *qp_colmaj_out_assign(ocp_qp_dims *qpdims, void *raw_memory) {
    return (qp_colmaj_out *) qp_sqp_out_assign(qpdims, raw_memory);
}

qp_colmaj_out *prepare_qp_colmaj_out(ocp_qp_dims *qpdims) {
    return (qp_colmaj_out *) prepare_qp_sqp_out(qpdims);
}

int sqp_in_calculate_size(sqp_dims *dims, ocp_qp_dims *qpdims) {
    int N = dims->N;
    //
    // Compute size of workspace
    //
    int size = sizeof(sqp_in);
    // x0
    size += dims->nx[0]*sizeof(double);
    // lb
    size += (N+1)*sizeof(double *);
    for (int i=0; i<=N; ++i) size += qpdims->nb[i]*sizeof(double);
    // ub
    size += (N+1)*sizeof(double *);
    for (int i=0; i<=N; ++i) size += qpdims->nb[i]*sizeof(double);
    // lg
    size += (N+1)*sizeof(double *);
    for (int i=0; i<=N; ++i) size += qpdims->ng[i]*sizeof(double);
    // ug
    size += (N+1)*sizeof(double *);
    for (int i=0; i<=N; ++i) size += qpdims->ng[i]*sizeof(double);
    // idxb
    size += (N+1)*sizeof(int *);
    for (int i=0; i<=N; ++i) size += qpdims->nb[i]*sizeof(int);

    return size;
}

sqp_in *sqp_in_assign(sqp_dims *dims, ocp_qp_dims *qpdims, void *raw_memory) {
    int N = dims->N;
    //
    // Assign workspace
    //
    char *c_ptr = raw_memory;

    sqp_in *out = (sqp_in *) c_ptr;
    c_ptr += sizeof(sqp_in);

    // x0
    out->x0 = (double *) c_ptr;
    c_ptr += dims->nx[0]*sizeof(double);
    // lb
    out->lb = (double **) c_ptr;
    c_ptr += (N+1)*sizeof(double *);
    for (int i=0; i<=N; ++i) {
        out->lb[i] = (double *) c_ptr;
        c_ptr += qpdims->nb[i]*sizeof(double);
    }
    // ub
    out->ub = (double **) c_ptr;
    c_ptr += (N+1)*sizeof(double *);
    for (int i=0; i<=N; ++i) {
        out->ub[i] = (double *) c_ptr;
        c_ptr += qpdims->nb[i]*sizeof(double);
    }
    // lg
    out->lg = (double **) c_ptr;
    c_ptr += (N+1)*sizeof(double *);
    for (int i=0; i<=N; ++i) {
        out->lg[i] = (double *) c_ptr;
        c_ptr += qpdims->ng[i]*sizeof(double);
    }
    // ug
    out->ug = (double **) c_ptr;
    c_ptr += (N+1)*sizeof(double *);
    for (int i=0; i<=N; ++i) {
        out->ug[i] = (double *) c_ptr;
        c_ptr += qpdims->ng[i]*sizeof(double);
    }
    // idxb
    out->idxb = (int **) c_ptr;
    c_ptr += (N+1)*sizeof(int *);
    for (int i=0; i<=N; ++i) {
        out->idxb[i] = (int *) c_ptr;
        c_ptr += qpdims->nb[i]*sizeof(int);
    }
    
    //
    // Check
    //
    assert(c_ptr == raw_memory + sqp_in_calculate_size(dims, qpdims) && "Assignment error.");

    return out;
}

sqp_in *prepare_sqp_in(sqp_dims *dims, ocp_qp_dims *qpdims) {
    //
    // Compute size of workspace
    //
    int size = sqp_in_calculate_size(dims, qpdims);

    //
    // Create workspace
    //
    char *raw_mem = (char *) malloc(size);

    //
    // Assign workspace
    //
    sqp_in *out = sqp_in_assign(dims, qpdims, raw_mem);

    return out;
}

int sqp_out_calculate_size(ocp_qp_dims *qpdims) {
    return qp_sqp_out_calculate_size(qpdims);
}

sqp_out *sqp_out_assign(ocp_qp_dims *qpdims, void *raw_memory) {
    return (sqp_out *) qp_sqp_out_assign(qpdims, raw_memory);
}

sqp_out *prepare_sqp_out(ocp_qp_dims *qpdims) {
    return (sqp_out *) prepare_qp_sqp_out(qpdims);
}

sqp_data *prepare_sqp_data(sqp_dims *dims, sqp_opts *opts, ocp_qp_dims *qpdims, qp_solver_config *qpconfig, void *qpopts) {
    int N = qpdims->N;
    //
    // Compute size of workspace
    //
    int size = sizeof(sqp_data);
    // rh_mem & lag_hess
    size += (N+1)*sizeof(regularize_hessian_mem *);
    size += (N+1)*sizeof(double *);
    for (int i=0; i<=N; ++i) {
        int n = dims->nx[i]+dims->nu[i];
        size += regularize_hessian_calculate_mem_size(n);
        size += n*n*sizeof(double);
    }
    // stage_cost
    size += (N+1)*sizeof(double);
    // colmaj_in
    size += qp_colmaj_in_calculate_size(qpdims);
    // colmaj_out
    size += qp_colmaj_out_calculate_size(qpdims);
    // colmaj_out_backup
    size += qp_colmaj_out_calculate_size(qpdims);
    // qpdims
    size += ocp_qp_dims_calculate_size(N);
    // qpin
    size += ocp_qp_in_calculate_size(qpconfig, qpdims);    
    // qpout
    size += ocp_qp_out_calculate_size(qpconfig, qpdims);
    // qpconfig
    size += sizeof(qp_solver_config);
    // qpopts
    // size += qpconfig->opts_calculate_size(qpconfig, qpdims);
    size += sizeof(void *);
    // qpmem
    size += qpconfig->memory_calculate_size(qpconfig, qpdims, qpopts);
    // qpwork
    size += qpconfig->workspace_calculate_size(qpconfig, qpdims, qpopts);
    // jk_memory
    size += (N+1)*sizeof(casadi_memory *);
    for (int i=0; i<N+1; ++i) {
        size += calculate_casadi_memory_size(opts->jk_fun[i], opts->jk_work[i], opts->jk_sparsity_in[i], opts->jk_sparsity_out[i]);
    }
    // gk_memory
    size += (N+1)*sizeof(casadi_memory *);
    for (int i=0; i<N+1; ++i) {
        size += calculate_casadi_memory_size(opts->gk_fun[i], opts->gk_work[i], opts->gk_sparsity_in[i], opts->gk_sparsity_out[i]);
    }
    // fk_memory
    size += N*sizeof(casadi_memory *);
    for (int i=0; i<N; ++i) {
        size += calculate_casadi_memory_size(opts->fk_fun[i], opts->fk_work[i], opts->fk_sparsity_in[i], opts->fk_sparsity_out[i]);
    }

    //
    // Create workspace
    //
    char *raw_mem = (char *) malloc(size);

    //
    // Assign workspace
    //
    char *c_ptr = raw_mem;

    sqp_data *out = (sqp_data *) c_ptr;
    c_ptr += sizeof(sqp_data);

    // rh_mem & lag_hess & interm_lag_hess
    out->rh_mem = (regularize_hessian_mem **) c_ptr;
    c_ptr += (N+1)*sizeof(regularize_hessian_mem *);
    out->lag_hess = (double **) c_ptr;
    c_ptr += (N+1)*sizeof(double *);
    for (int i=0; i<=N; ++i) {
        int n = dims->nx[i]+dims->nu[i];
        out->rh_mem[i] = regularize_hessian_assign_mem(n, c_ptr);
        c_ptr += regularize_hessian_calculate_mem_size(n);
        out->lag_hess[i] = (double *) c_ptr;
        c_ptr += n*n*sizeof(double);
    }
    // stage_cost
    out->stage_cost = (double *) c_ptr;
    c_ptr += (N+1)*sizeof(double);
    // colmaj_in
    out->colmaj_in = qp_colmaj_in_assign(qpdims, c_ptr);
    c_ptr += qp_colmaj_in_calculate_size(qpdims);
    // colmaj_out
    out->colmaj_out = qp_colmaj_out_assign(qpdims, c_ptr);
    c_ptr += qp_colmaj_out_calculate_size(qpdims);
    // colmaj_out_backup
    out->colmaj_out_backup = qp_colmaj_out_assign(qpdims, c_ptr);
    c_ptr += qp_colmaj_out_calculate_size(qpdims);
    // qpdims
    out->qpdims = ocp_qp_dims_assign(N, c_ptr);
    c_ptr += ocp_qp_dims_calculate_size(N);
    // qpin
    out->qpin = ocp_qp_in_assign(qpconfig, qpdims, c_ptr);
    c_ptr += ocp_qp_in_calculate_size(qpconfig, qpdims);    
    // qpout
    out->qpout = ocp_qp_out_assign(qpconfig, qpdims, c_ptr);
    c_ptr += ocp_qp_out_calculate_size(qpconfig, qpdims);
    // qpconfig
    out->qpconfig = (qp_solver_config *) c_ptr;
    c_ptr += sizeof(qp_solver_config);
    // qpopts: DON'T LIKE THIS, BUT APPARENTLY THIS IS ACADOS STYLE
    // out->qpopts = qpconfig->opts_assign(qpconfig, qpdims, c_ptr);
    // c_ptr += qpconfig->opts_calculate_size(qpconfig, qpdims);
    out->qpopts = qpopts;
    c_ptr += sizeof(void *);
    // qpmem
    out->qpmem = qpconfig->memory_assign(qpconfig, qpdims, qpopts, c_ptr);
    c_ptr += qpconfig->memory_calculate_size(qpconfig, qpdims, qpopts);
    // qpwork
    out->qpwork = (void *) c_ptr;
    c_ptr += qpconfig->workspace_calculate_size(qpconfig, qpdims, qpopts);
    // jk_memory & gk_memory
    out->jk_memory = (casadi_memory **) c_ptr;
    c_ptr += (N+1)*sizeof(casadi_memory *);
    out->gk_memory = (casadi_memory **) c_ptr;
    c_ptr += (N+1)*sizeof(casadi_memory *);
    for (int i=0; i<N+1;++i) {
        out->jk_memory[i] = assign_casadi_memory(opts->jk_fun[i], opts->jk_work[i], opts->jk_sparsity_in[i], opts->jk_sparsity_out[i], c_ptr);
        c_ptr += calculate_casadi_memory_size(opts->jk_fun[i], opts->jk_work[i], opts->jk_sparsity_in[i], opts->jk_sparsity_out[i]);
        out->gk_memory[i] = assign_casadi_memory(opts->gk_fun[i], opts->gk_work[i], opts->gk_sparsity_in[i], opts->gk_sparsity_out[i], c_ptr);
        c_ptr += calculate_casadi_memory_size(opts->gk_fun[i], opts->gk_work[i], opts->gk_sparsity_in[i], opts->gk_sparsity_out[i]);
    }
    // fk_memory
    out->fk_memory = (casadi_memory **) c_ptr;
    c_ptr += N*sizeof(casadi_memory *);
    for (int i=0; i<N;++i) {
        out->fk_memory[i] = assign_casadi_memory(opts->fk_fun[i], opts->fk_work[i], opts->fk_sparsity_in[i], opts->fk_sparsity_out[i], c_ptr);
        c_ptr += calculate_casadi_memory_size(opts->fk_fun[i], opts->fk_work[i], opts->fk_sparsity_in[i], opts->fk_sparsity_out[i]);
    }

    //
    // Check
    //
    assert(c_ptr == raw_mem + size && "Assignment error.");

    //
    // Configure
    //
    // qpdims
    out->qpdims->N = qpdims->N;
    for (int i=0; i<=N; ++i) {
        out->qpdims->nx[i] = qpdims->nx[i];
        out->qpdims->nu[i] = qpdims->nu[i];
        out->qpdims->nbx[i] = qpdims->nbx[i];
        out->qpdims->nbu[i] = qpdims->nbu[i];
        out->qpdims->nb[i] = qpdims->nb[i];
        out->qpdims->ng[i] = qpdims->ng[i];
        out->qpdims->ns[i] = qpdims->ns[i];
    }
    // qpconfig
    *out->qpconfig = *qpconfig;

    return out;
}

int sqp(sqp_dims *dims, sqp_in *in, sqp_out *out, sqp_opts *opts, sqp_data *data) {
    int sqp_its = opts->its;
    int status = 0;

    ocp_qp_dims *qpdims = data->qpdims;
    int N = qpdims->N;


    // Bounds
    for (int ii=0; ii<=N; ++ii) {
        memcpy(data->colmaj_in->idxb[ii], in->idxb[ii], qpdims->nb[ii]*sizeof(int));
        memcpy(data->colmaj_in->lb[ii], in->lb[ii], qpdims->nb[ii]*sizeof(double));
        memcpy(data->colmaj_in->ub[ii], in->ub[ii], qpdims->nb[ii]*sizeof(double));
        for (int iii=0; iii<qpdims->nb[ii]; ++iii) {
            int idxbi = data->colmaj_in->idxb[ii][iii];
            if (idxbi < qpdims->nu[ii]) {
                data->colmaj_in->lb[ii][iii] -= out->u[ii][idxbi];
                data->colmaj_in->ub[ii][iii] -= out->u[ii][idxbi];
            } else {
                idxbi -= qpdims->nu[ii];
                data->colmaj_in->lb[ii][iii] -= out->x[ii][idxbi];
                data->colmaj_in->ub[ii][iii] -= out->x[ii][idxbi];
            }
        }
    }



    for (int i=0; i<sqp_its; ++i) {
        // Reset Hessian
        for (int ii=0; ii<=N; ++ii) {
            int n = dims->nx[ii]+dims->nu[ii];
            memset(data->lag_hess[ii], 0, n*n*sizeof(double));
        }
        


        // Evaluate dynamics
        //  first stage
        eval_fk(in->x0, out->u[0], out->pi[0], data->lag_hess[0], in->p, data->lag_hess[0], data->colmaj_in->A[0], data->colmaj_in->B[0], data->colmaj_in->b[0], data->fk_memory[0]);
        for (int iii=0; iii<qpdims->nx[1]; ++iii) data->colmaj_in->b[0][iii] -= out->x[1][iii];
        //  middle stages
        for (int ii=1; ii<N; ++ii) {
            eval_fk(out->x[ii], out->u[ii], out->pi[ii], data->lag_hess[ii], in->p, data->lag_hess[ii], data->colmaj_in->A[ii], data->colmaj_in->B[ii], data->colmaj_in->b[ii], data->fk_memory[ii]);
            for (int iii=0; iii<qpdims->nx[ii+1]; ++iii) data->colmaj_in->b[ii][iii] -= out->x[ii+1][iii];
        }



        // Evaluate nonlinear constraints
        //  first stage
        eval_gk(in->x0, out->u[0], out->lam_lg[0], out->lam_ug[0], data->lag_hess[0], in->p, data->lag_hess[0], data->colmaj_in->C[0], data->colmaj_in->D[0], data->colmaj_in->lg[0], data->colmaj_in->ug[0], data->gk_memory[0]); //LOWER AND UPPER lagrange multipliers.
        for (int iii=0; iii<qpdims->ng[0]; ++iii) {
            data->colmaj_in->lg[0][iii] += in->lg[0][iii];
            data->colmaj_in->ug[0][iii] += in->ug[0][iii];
        }
        //  middle stages
        for (int ii=1; ii<=N; ++ii) {
            eval_gk(out->x[ii], out->u[ii], out->lam_lg[ii], out->lam_ug[ii], data->lag_hess[ii], in->p, data->lag_hess[ii], data->colmaj_in->C[ii], data->colmaj_in->D[ii], data->colmaj_in->lg[ii], data->colmaj_in->ug[ii], data->gk_memory[ii]); //LOWER AND UPPER lagrange multipliers.
            for (int iii=0; iii<qpdims->ng[ii]; ++iii) {
                data->colmaj_in->lg[ii][iii] += in->lg[ii][iii];
                data->colmaj_in->ug[ii][iii] += in->ug[ii][iii];
            }
        }



        // Evaluate cost
        //  first stage
        eval_jk(in->x0, out->u[0], data->lag_hess[0], in->p, data->lag_hess[0], data->colmaj_in->q[0], data->colmaj_in->r[0], data->stage_cost, data->jk_memory[0]);
        //  middle stages
        for (int ii=1; ii<=N; ++ii) {
            eval_jk(out->x[ii], out->u[ii], data->lag_hess[ii], in->p, data->lag_hess[ii], data->colmaj_in->q[ii], data->colmaj_in->r[ii], data->stage_cost + ii, data->jk_memory[ii]);
        }



        // Compute Hessian of Lagrangian and regularize
        for (int ii=0; ii<=N; ++ii) {
            int n = dims->nx[ii]+dims->nu[ii];
            // Regularize
            if (opts->eh) {
                regularize_hessian(data->lag_hess[ii], n, opts->reg_eps, data->rh_mem[ii]);
            }

            // Write to QP def
            for (int iii=0; iii<qpdims->nu[ii]; ++iii) {
                for (int iiii=0; iiii<qpdims->nx[ii]; ++iiii) data->colmaj_in->S[ii][iiii*qpdims->nu[ii]+iii] = data->lag_hess[ii][iiii*n+dims->nx[ii]+iii];
                for (int iiii=0; iiii<qpdims->nu[ii]; ++iiii) data->colmaj_in->R[ii][iiii*qpdims->nu[ii]+iii] = data->lag_hess[ii][(iiii+dims->nx[ii])*n+dims->nx[ii]+iii];
            }
            for (int iii=0; iii<qpdims->nx[ii]; ++iii) {
                for (int iiii=0; iiii<qpdims->nx[ii]; ++iiii) data->colmaj_in->Q[ii][iiii*qpdims->nx[ii]+iii] = data->lag_hess[ii][iiii*n+iii];
            }
        }



        // Pack QP
        d_cvt_colmaj_to_ocp_qp(data->colmaj_in->A, data->colmaj_in->B, data->colmaj_in->b, data->colmaj_in->Q, data->colmaj_in->S, data->colmaj_in->R, data->colmaj_in->q, data->colmaj_in->r, data->colmaj_in->idxb, data->colmaj_in->lb, data->colmaj_in->ub, data->colmaj_in->C, data->colmaj_in->D, data->colmaj_in->lg, data->colmaj_in->ug, NULL, NULL, NULL, NULL, NULL, NULL, NULL, data->qpin);



        // Solve QP
        status = data->qpconfig->evaluate(data->qpconfig, data->qpin, data->qpout, data->qpopts, data->qpmem, data->qpwork);

        if (status!=0) {
            if (i > 0) {
                printf("Non-zero QP-solver status in iteration %d of the SQP method: %d --  Backtracking...\n", i, status);
                
                
                
                // Backtrack solution
                double backtrack = 0.1;
                for (int ii=0;ii<=N;++ii) {
                    for (int iii=0;iii<qpdims->nx[ii];++iii) {
                        out->x[ii][iii] = data->colmaj_out_backup->x[ii][iii] + backtrack*data->colmaj_out->x[ii][iii];
                        data->colmaj_out->x[ii][iii] = backtrack*data->colmaj_out->x[ii][iii];
                    }
                    for (int iii=0;iii<qpdims->nu[ii];++iii) {
                        out->u[ii][iii] = data->colmaj_out_backup->u[ii][iii] + backtrack*data->colmaj_out->u[ii][iii];
                        data->colmaj_out->u[ii][iii] = backtrack*data->colmaj_out->u[ii][iii];
                    }
                    for (int iii=0;iii<qpdims->nb[ii];++iii) out->lam_lb[ii][iii] = data->colmaj_out_backup->lam_lb[ii][iii];
                    for (int iii=0;iii<qpdims->nb[ii];++iii) out->lam_ub[ii][iii] = data->colmaj_out_backup->lam_ub[ii][iii];
                    for (int iii=0;iii<qpdims->ng[ii];++iii) out->lam_lg[ii][iii] = data->colmaj_out_backup->lam_lg[ii][iii];
                    for (int iii=0;iii<qpdims->ng[ii];++iii) out->lam_ug[ii][iii] = data->colmaj_out_backup->lam_ug[ii][iii];
                }
                for (int ii=0;ii<N;++ii) {
                    for (int iii=0;iii<qpdims->nx[ii+1];++iii) out->pi[ii][iii] = data->colmaj_out_backup->pi[ii][iii];
                }
            } else {
                printf("Non-zero QP-solver status in iteration %d of the SQP method: %d --  In first iteration (what can I do?!)...\n", i, status);

                // Unpack QP
                d_cvt_ocp_qp_sol_to_colmaj(data->qpout, data->colmaj_out->u, data->colmaj_out->x, NULL, NULL, data->colmaj_out->pi, data->colmaj_out->lam_lb, data->colmaj_out->lam_ub, data->colmaj_out->lam_lg, data->colmaj_out->lam_ug, NULL, NULL);



                // Increment solution
                for (int ii=0;ii<=N;++ii) {
                    for (int iii=0;iii<qpdims->nx[ii];++iii) data->colmaj_out_backup->x[ii][iii] = out->x[ii][iii];
                    for (int iii=0;iii<qpdims->nx[ii];++iii) out->x[ii][iii] += data->colmaj_out->x[ii][iii];
                    for (int iii=0;iii<qpdims->nu[ii];++iii) data->colmaj_out_backup->u[ii][iii] = out->u[ii][iii];
                    for (int iii=0;iii<qpdims->nu[ii];++iii) out->u[ii][iii] += data->colmaj_out->u[ii][iii];
                    for (int iii=0;iii<qpdims->nb[ii];++iii) data->colmaj_out_backup->lam_lb[ii][iii] = out->lam_lb[ii][iii];
                    for (int iii=0;iii<qpdims->nb[ii];++iii) out->lam_lb[ii][iii] = data->colmaj_out->lam_lb[ii][iii];
                    for (int iii=0;iii<qpdims->nb[ii];++iii) data->colmaj_out_backup->lam_ub[ii][iii] = out->lam_ub[ii][iii];
                    for (int iii=0;iii<qpdims->nb[ii];++iii) out->lam_ub[ii][iii] = data->colmaj_out->lam_ub[ii][iii];
                    for (int iii=0;iii<qpdims->ng[ii];++iii) data->colmaj_out_backup->lam_lg[ii][iii] = out->lam_lg[ii][iii];
                    for (int iii=0;iii<qpdims->ng[ii];++iii) out->lam_lg[ii][iii] = data->colmaj_out->lam_lg[ii][iii];
                    for (int iii=0;iii<qpdims->ng[ii];++iii) data->colmaj_out_backup->lam_ug[ii][iii] = out->lam_ug[ii][iii];
                    for (int iii=0;iii<qpdims->ng[ii];++iii) out->lam_ug[ii][iii] = data->colmaj_out->lam_ug[ii][iii];

#if (CAPCONTROLS==1)
                    // Cap controls
                    for (int iii=0;iii<qpdims->nb[ii];++iii) {
                        int idxbi = in->idxb[ii][iii];
                        if (idxbi < qpdims->nu[ii]) {
                            out->u[ii][idxbi] = fmin(in->ub[ii][iii], fmax(in->lb[ii][iii], out->u[ii][idxbi]));
                        }
                    }
#endif
                }
                for (int ii=0;ii<N;++ii) {
                    for (int iii=0;iii<qpdims->nx[ii+1];++iii) data->colmaj_out_backup->pi[ii][iii] = out->pi[ii][iii];
                    for (int iii=0;iii<qpdims->nx[ii+1];++iii) out->pi[ii][iii] = data->colmaj_out->pi[ii][iii];
                }
            }
        } else {
            // Unpack QP
            d_cvt_ocp_qp_sol_to_colmaj(data->qpout, data->colmaj_out->u, data->colmaj_out->x, NULL, NULL, data->colmaj_out->pi, data->colmaj_out->lam_lb, data->colmaj_out->lam_ub, data->colmaj_out->lam_lg, data->colmaj_out->lam_ug, NULL, NULL);



            // Increment solution
            for (int ii=0;ii<=N;++ii) {
                for (int iii=0;iii<qpdims->nx[ii];++iii) data->colmaj_out_backup->x[ii][iii] = out->x[ii][iii];
                for (int iii=0;iii<qpdims->nx[ii];++iii) out->x[ii][iii] += data->colmaj_out->x[ii][iii];
                for (int iii=0;iii<qpdims->nu[ii];++iii) data->colmaj_out_backup->u[ii][iii] = out->u[ii][iii];
                for (int iii=0;iii<qpdims->nu[ii];++iii) out->u[ii][iii] += data->colmaj_out->u[ii][iii];
                for (int iii=0;iii<qpdims->nb[ii];++iii) data->colmaj_out_backup->lam_lb[ii][iii] = out->lam_lb[ii][iii];
                for (int iii=0;iii<qpdims->nb[ii];++iii) out->lam_lb[ii][iii] = data->colmaj_out->lam_lb[ii][iii];
                for (int iii=0;iii<qpdims->nb[ii];++iii) data->colmaj_out_backup->lam_ub[ii][iii] = out->lam_ub[ii][iii];
                for (int iii=0;iii<qpdims->nb[ii];++iii) out->lam_ub[ii][iii] = data->colmaj_out->lam_ub[ii][iii];
                for (int iii=0;iii<qpdims->ng[ii];++iii) data->colmaj_out_backup->lam_lg[ii][iii] = out->lam_lg[ii][iii];
                for (int iii=0;iii<qpdims->ng[ii];++iii) out->lam_lg[ii][iii] = data->colmaj_out->lam_lg[ii][iii];
                for (int iii=0;iii<qpdims->ng[ii];++iii) data->colmaj_out_backup->lam_ug[ii][iii] = out->lam_ug[ii][iii];
                for (int iii=0;iii<qpdims->ng[ii];++iii) out->lam_ug[ii][iii] = data->colmaj_out->lam_ug[ii][iii];

#if (CAPCONTROLS==1)
                // Cap controls
                for (int iii=0;iii<qpdims->nb[ii];++iii) {
                    int idxbi = in->idxb[ii][iii];
                    if (idxbi < qpdims->nu[ii]) {
                        out->u[ii][idxbi] = fmin(in->ub[ii][iii], fmax(in->lb[ii][iii], out->u[ii][idxbi]));
                    }
                }
#endif
            }
            for (int ii=0;ii<N;++ii) {
                for (int iii=0;iii<qpdims->nx[ii+1];++iii) data->colmaj_out_backup->pi[ii][iii] = out->pi[ii][iii];
                for (int iii=0;iii<qpdims->nx[ii+1];++iii) out->pi[ii][iii] = data->colmaj_out->pi[ii][iii];
            }
        }

        double kkt = compute_kkt(data);

        if (kkt < opts->tol) {
            break;
        }

        // printf("KKT value: %.5e\n", compute_kkt(data));
    }

    return status;
}

double compute_kkt(sqp_data *data) {
    double kkt = 0;

    ocp_qp_dims *qpdims = data->qpdims;

    for(int i=0; i<=qpdims->N; ++i) {
        for(int ii=0; ii<qpdims->nx[i]; ++ii) kkt += fabs(data->colmaj_in->q[i][ii]*data->colmaj_out->x[i][ii]);
        for(int ii=0; ii<qpdims->nu[i]; ++ii) kkt += fabs(data->colmaj_in->r[i][ii]*data->colmaj_out->u[i][ii]);
    }
    
    for(int i=0; i<qpdims->N; ++i) {
        for(int ii=0; ii<qpdims->nx[i+1]; ++ii) kkt += fabs(data->colmaj_in->b[i][ii]*data->colmaj_out->pi[i][ii]);
    }

    for(int i=0; i<=qpdims->N; ++i) {
        for(int ii=0; ii<qpdims->ng[i]; ++ii) kkt += fabs(-data->colmaj_in->lg[i][ii]*data->colmaj_out->lam_lg[i][ii]);
        for(int ii=0; ii<qpdims->ng[i]; ++ii) kkt += fabs(data->colmaj_in->ug[i][ii]*data->colmaj_out->lam_ug[i][ii]);
    }

    return kkt;
}
