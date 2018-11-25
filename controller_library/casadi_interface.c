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

#include "casadi_interface.h"

#include <assert.h>
#include <stdlib.h>

/////////////////////////////////////////
//        FORWARD DECLARATIONS         //
/////////////////////////////////////////

int nnz(const int *sparsity);
void densify(const double *sparse_in, double *dense_out, const int *sparsity);
void sparsify(const double *dense_in, double *sparse_out, const int *sparsity);

/////////////////////////////////////////
//         CONVENIENCE WRAPPERS        //
/////////////////////////////////////////

int eval_jk(const double *_Xk, const double *_Uk, const double *_HESSk, const double *_P, double *HESSk_, double *qk_, double *rk_, double *jk_, casadi_memory *mem) {
    int status;

    //
    // Write to arg
    //
    sparsify(_Xk, (double *)mem->arg[0], mem->sparsity_in(0));
    sparsify(_Uk, (double *)mem->arg[1], mem->sparsity_in(1));
    sparsify(_HESSk, (double *)mem->arg[2], mem->sparsity_in(2));
    sparsify(_P, (double *)mem->arg[3], mem->sparsity_in(3));

    //
    // Call routine
    //
    status = mem->fun(mem->arg, mem->res, mem->iw, mem->w, 0);

    //
    // Write from res
    //
    densify(mem->res[0], HESSk_, mem->sparsity_out(0));
    densify(mem->res[1], qk_, mem->sparsity_out(1));
    densify(mem->res[2], rk_, mem->sparsity_out(2));
    densify(mem->res[3], jk_, mem->sparsity_out(3));

    return status;
}

int eval_gk(const double *_Xk, const double *_Uk, const double *_LMUk, const double *_UMUk, const double *_HESSk, const double *_P, double *HESSk_, double *Cxk_, double *Cuk_, double *lck_, double *uck_, casadi_memory *mem) {
    int status;

    //
    // Write to arg
    //
    sparsify(_Xk, (double *)mem->arg[0], mem->sparsity_in(0));
    sparsify(_Uk, (double *)mem->arg[1], mem->sparsity_in(1));
    sparsify(_LMUk, (double *)mem->arg[2], mem->sparsity_in(2));
    sparsify(_UMUk, (double *)mem->arg[3], mem->sparsity_in(3));
    sparsify(_HESSk, (double *)mem->arg[4], mem->sparsity_in(4));
    sparsify(_P, (double *)mem->arg[5], mem->sparsity_in(5));

    //
    // Call routine
    //
    status = mem->fun(mem->arg, mem->res, mem->iw, mem->w, 0);

    //
    // Write from res
    //
    densify(mem->res[0], HESSk_, mem->sparsity_out(0));
    densify(mem->res[1], Cxk_, mem->sparsity_out(1));
    densify(mem->res[2], Cuk_, mem->sparsity_out(2));
    densify(mem->res[3], lck_, mem->sparsity_out(3));
    densify(mem->res[4], uck_, mem->sparsity_out(4));

    return status;
}

int eval_fk(const double *_Xk, const double *_Uk, const double *_PIk, const double *_HESSk, const double *_P, double *HESSk_, double *A_, double *B_, double *b_, casadi_memory *mem) {
    int status;

    //
    // Write to arg
    //
    sparsify(_Xk, (double *)mem->arg[0], mem->sparsity_in(0));
    sparsify(_Uk, (double *)mem->arg[1], mem->sparsity_in(1));
    sparsify(_PIk, (double *)mem->arg[2], mem->sparsity_in(2));
    sparsify(_HESSk, (double *)mem->arg[3], mem->sparsity_in(3));
    sparsify(_P, (double *)mem->arg[4], mem->sparsity_in(4));

    //
    // Call routine
    //
    status = mem->fun(mem->arg, mem->res, mem->iw, mem->w, 0);

    //
    // Write from res
    //
    densify(mem->res[0], HESSk_, mem->sparsity_out(0));
    densify(mem->res[1], A_, mem->sparsity_out(1));
    densify(mem->res[2], B_, mem->sparsity_out(2));
    densify(mem->res[3], b_, mem->sparsity_out(3));

    return status;
}

/////////////////////////////////////////
//           HELPER ROUTINES           //
/////////////////////////////////////////

int calculate_casadi_memory_size(int (*cas_fun)(), int (*cas_work)(), const int *(*cas_sparsity_in)(), const int *(*cas_sparsity_out)()) {
    //
    // Compute size of casadi_memory
    //
    int bytes = 0;

    int sz_arg, sz_res, sz_iw, sz_w;
    cas_work(&sz_arg,&sz_res,&sz_iw,&sz_w);

    // casadi_memory-struct
    bytes += sizeof(casadi_memory);
    // arg
    bytes += sz_arg*sizeof(const double *);
    for (int i=0;i<sz_arg;++i) {
        const int *sparsity = cas_sparsity_in(i);
        bytes += nnz(sparsity)*sizeof(const double);
    }
    // res
    bytes += sz_res*sizeof(double *);
    for (int i=0;i<sz_res;++i) {
        const int *sparsity = cas_sparsity_out(i);
        bytes += nnz(sparsity)*sizeof(double);
    }
    // iw
    bytes += sz_iw*sizeof(int);
    // w
    bytes += sz_w*sizeof(double);

    return bytes;
}

casadi_memory *assign_casadi_memory(int (*cas_fun)(), int (*cas_work)(), const int *(*cas_sparsity_in)(), const int *(*cas_sparsity_out)(), void* raw_memory) {
    int sz_arg, sz_res, sz_iw, sz_w;
    cas_work(&sz_arg, &sz_res, &sz_iw, &sz_w);

    //
    // Assign casadi_memory 
    //
    char *c_ptr = (char *) raw_memory;

    // casadi_memory-struct
    casadi_memory *out = (casadi_memory *) c_ptr;
    c_ptr += sizeof(casadi_memory);
    // arg
    out->arg = (const double **) c_ptr;
    c_ptr += sz_arg*sizeof(const double *);
    for (int i=0;i<sz_arg;++i) {
        const int *sparsity = cas_sparsity_in(i);
        out->arg[i] = (const double *) c_ptr;
        c_ptr += nnz(sparsity)*sizeof(const double);
    }
    // res
    out->res = (double **) c_ptr;
    c_ptr += sz_res*sizeof(double *);
    for (int i=0;i<sz_res;++i) {
        const int *sparsity = cas_sparsity_out(i);
        out->res[i] = (double *) c_ptr;
        c_ptr += nnz(sparsity)*sizeof(double);
    }
    // iw
    out->iw = (int *) c_ptr;
    c_ptr += sz_iw*sizeof(int);
    // w
    out->w = (double *) c_ptr;
    c_ptr += sz_w*sizeof(double);

    assert(c_ptr == raw_memory + calculate_casadi_memory_size(cas_fun, cas_work, cas_sparsity_in, cas_sparsity_out) && "Memory assignment mismatch.");

    //
    // Configure
    //
    out->work = cas_work;
    out->sparsity_in = cas_sparsity_in;
    out->sparsity_out = cas_sparsity_out;
    out->fun = cas_fun;

    return out;
}

casadi_memory *create_casadi_memory(int (*cas_fun)(), int (*cas_work)(), const int *(*cas_sparsity_in)(), const int *(*cas_sparsity_out)()) {
    //
    // Compute size of casadi_memory
    //
    int bytes = calculate_casadi_memory_size(cas_fun, cas_work, cas_sparsity_in, cas_sparsity_out);

    //
    // Create casadi_memory
    //
    void *raw_mem = malloc(bytes);

    //
    // Assign casadi_memory 
    //
    casadi_memory *out = assign_casadi_memory(cas_fun, cas_work, cas_sparsity_in, cas_sparsity_out, raw_mem);

    return out;
}

int nnz(const int *sparsity)
{
    int nnz = 0;
    if (sparsity != NULL) {
        const int nrow = sparsity[0];
        const int ncol = sparsity[1];
        const int dense = sparsity[2];
        if (dense) {
            nnz = nrow * ncol;
        } else {
            const int *colind = sparsity + 2;
            for (int i = 0; i < ncol; ++i) {
                nnz += colind[i + 1] - colind[i];
            }
        }
    }
    return nnz;
}

void densify(const double *sparse_in, double *dense_out, const int *sparsity)
{
    if (dense_out != NULL) {
        const int nrow = sparsity[0];
        const int ncol = sparsity[1];
        const int dense = sparsity[2];
        if (dense) {
            for (int i = 0; i < ncol * nrow; i++) dense_out[i] = sparse_in[i];
        } else {
            // Fill with zeros
            for (int i = 0; i < ncol; i++) {
                for (int j = 0; j < nrow; j++) dense_out[i * nrow + j] = 0;
            }
            // Additional data
            const double *x = sparse_in;
            const int *colind = sparsity + 2, *row = sparsity + ncol + 3;
            // Copy nonzeros
            for (int i = 0; i < ncol; ++i) {
                for (int el = colind[i]; el != colind[i + 1]; ++el) {
                    dense_out[row[el] + i * nrow] = *x++;
                }
            }
        }
    }
}

void sparsify(const double *dense_in, double *sparse_out, const int *sparsity)
{
    const int nrow = sparsity[0];
    const int ncol = sparsity[1];
    const int dense = sparsity[2];
    if (dense) {
        for (int i = 0; i < ncol * nrow; i++) sparse_out[i] = dense_in[i];
    } else {
        // Additional data
        double *x = sparse_out;
        const int *colind = sparsity + 2, *row = sparsity + ncol + 3;
        // Copy nonzeros
        for (int i = 0; i < ncol; ++i) {
            for (int el = colind[i]; el != colind[i + 1]; ++el) {
                *x++ = dense_in[row[el] + i * nrow];
            }
        }
    }
}
