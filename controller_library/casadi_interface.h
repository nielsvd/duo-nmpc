#ifndef CASADI_INTERFACE_H_
#define CASADI_INTERFACE_H_

typedef struct casadi_memory_ {
    int *iw;
    double *w;
    const double **arg;
    double **res;
    int (*work)(int *, int *, int *, int *);
    const int *(*sparsity_in)(int);
    const int *(*sparsity_out)(int);
    int (*fun)(const double **, double **, int *, double *, void *);
} casadi_memory;

// Convenience wrappers
int eval_jk(const double *_Xk, const double *_Uk, const double *_HESSk, const double *_P, double *HESSk_, double *qk_, double *rk_, double *jk_, casadi_memory *mem);
int eval_gk(const double *_Xk, const double *_Uk, const double *_LMUk, const double *_UMUk, const double *_HESSk, const double *_P, double *HESSk_, double *Cxk_, double *Cuk_, double *lck_, double *uck_, casadi_memory *mem);
int eval_fk(const double *_Xk, const double *_Uk, const double *_PIk, const double *_HESSk, const double *_P, double *HESSk_, double *A_, double *B_, double *b_, casadi_memory *mem);

// Helper routines
int calculate_casadi_memory_size(int (*cas_fun)(), int (*cas_work)(), const int *(*cas_sparsity_in)(), const int *(*cas_sparsity_out)());
casadi_memory *assign_casadi_memory(int (*cas_fun)(), int (*cas_work)(), const int *(*cas_sparsity_in)(), const int *(*cas_sparsity_out)(), void *raw_memory);
casadi_memory *create_casadi_memory(int (*cas_fun)(), int (*cas_work)(), const int *(*cas_sparsity_in)(), const int *(*cas_sparsity_out)());

#endif // CASADI_INTERFACE_H_