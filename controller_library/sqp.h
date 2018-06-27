#ifndef SQP_H_
#define SQP_H_


#include <acados/ocp_qp/ocp_qp_common.h>
#include <acados_c/ocp_qp_interface.h>
#include <acados_c/options_interface.h>

#include "regularize_hessian.h"
#include "casadi_interface.h"

typedef struct sqp_dims_ {
    int N;
    int *nx;
    int *nu;
    int *ng;
    int *nb;
} sqp_dims;

typedef struct sqp_opts_ {
    int its;
    int eh;
    double tol;
    double reg_eps;
    // casadi functions
    int (**jk_work)(int*, int*, int*, int*);
    const int *(**jk_sparsity_in)(int);
    const int *(**jk_sparsity_out)(int);
    int (**jk_fun)(const double**, double**, int*, double*, void*);
    int (**gk_work)(int*, int*, int*, int*);
    const int *(**gk_sparsity_in)(int);
    const int *(**gk_sparsity_out)(int);
    int (**gk_fun)(const double**, double**, int*, double*, void*);
    int (**fk_work)(int*, int*, int*, int*);
    const int *(**fk_sparsity_in)(int);
    const int *(**fk_sparsity_out)(int);
    int (**fk_fun)(const double**, double**, int*, double*, void*);
} sqp_opts;

typedef struct qp_colmaj_in_ {
    double **A;
    double **B;
    double **b;
    double **Q;
    double **S;
    double **R;
    double **q;
    double **r;
    int **idxb;
    double **lb;
    double **ub;
    double **C;
    double **D;
    double **lg;
    double **ug;
} qp_colmaj_in;

typedef struct qp_colmaj_out_ {
    double **x;
    double **u;
    double **lam_lb;
    double **lam_ub;
    double **lam_lg;
    double **lam_ug;
    double **pi;
} qp_colmaj_out;

typedef struct sqp_in_ {
    const double *x0;
    const double *p;
    int **idxb;
    double **lb;
    double **ub;
    double **lg;
    double **ug;
} sqp_in;

typedef struct sqp_out_ {
    double **x;
    double **u;
    double **lam_lb;
    double **lam_ub;
    double **lam_lg;
    double **lam_ug;
    double **pi;
} sqp_out;

typedef struct sqp_data_ {
    regularize_hessian_mem **rh_mem;
    double **lag_hess;
    double *stage_cost;
    qp_colmaj_in *colmaj_in;
    qp_colmaj_out *colmaj_out;
    qp_colmaj_out *colmaj_out_backup;
    ocp_qp_dims *qpdims;
    ocp_qp_in *qpin;
    ocp_qp_out *qpout;
    qp_solver_config *qpconfig;
    void *qpopts;
    void *qpmem;
    void *qpwork;
    casadi_memory **jk_memory;
    casadi_memory **gk_memory;
    casadi_memory **fk_memory;
} sqp_data;

sqp_dims *prepare_sqp_dims(int N);
sqp_opts *prepare_sqp_opts(sqp_dims *dims);
int qp_colmaj_in_calculate_size(ocp_qp_dims *qpdims);
qp_colmaj_in *qp_colmaj_in_assign(ocp_qp_dims *qpdims, void *raw_memory);
qp_colmaj_in *prepare_qp_colmaj_in(ocp_qp_dims *qpdims);
int qp_colmaj_out_calculate_size(ocp_qp_dims *qpdims);
qp_colmaj_out *qp_colmaj_out_assign(ocp_qp_dims *qpdims, void *raw_memory);
qp_colmaj_out *prepare_qp_colmaj_out(ocp_qp_dims *qpdims);
int sqp_in_calculate_size(sqp_dims *dims, ocp_qp_dims *qpdims);
sqp_in *sqp_in_assign(sqp_dims *dims, ocp_qp_dims *qpdims, void *raw_memory);
sqp_in *prepare_sqp_in(sqp_dims *dims, ocp_qp_dims *qpdims);
int sqp_out_calculate_size(ocp_qp_dims *qpdims);
sqp_out *sqp_out_assign(ocp_qp_dims *qpdims, void *raw_memory);
sqp_out *prepare_sqp_out(ocp_qp_dims *qpdims);
sqp_data *prepare_sqp_data(sqp_dims *dims, sqp_opts *opts, ocp_qp_dims *qpdims, qp_solver_config *qpconfig, void *qpopts);
int sqp(sqp_dims *dims, sqp_in *in, sqp_out *out, sqp_opts *opts, sqp_data *data);
double compute_kkt(sqp_data *data);

#endif // SQP_H_