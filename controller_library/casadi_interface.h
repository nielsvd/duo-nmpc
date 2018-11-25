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