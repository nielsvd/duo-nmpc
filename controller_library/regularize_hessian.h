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

#ifndef REGULARIZE_HESSIAN_H_
#define REGULARIZE_HESSIAN_H_

typedef struct regularize_hessian_mem_ {
    double *Q;  // Eigenvectors
    double *Ld; // Eigenvalues
    double *QL;
    int (*blas_daxpy)(int *n, double *da, double *dx, int *incx, double *dy, int *incy);
    void (*blas_dsyev)(const char *jobz, const char *uplo, int *n, double *a, int *lda, double *w, double *work, int *lwork, int *info);
    int (*blas_dgemm)(const char *transa, const char *transb, int *m, int *n, int *k, double *alpha, double *a, int *lda, double *b, int *ldb, double *beta, double *c, int *ldc);
} regularize_hessian_mem;

int regularize_hessian_calculate_mem_size(int n);
regularize_hessian_mem *regularize_hessian_assign_mem(int n, void *raw_memory);
regularize_hessian_mem *regularize_hessian_create_mem(int n);
int regularize_hessian(double *H, int n, double eps, regularize_hessian_mem *mem);

#endif // REGULARIZE_HESSIAN_H_