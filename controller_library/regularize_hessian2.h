#ifndef REGULARIZE_HESSIAN2_H_
#define REGULARIZE_HESSIAN2_H_

typedef struct regularize_hessian2_mem_ {
	double *V; // Eigenvectors
	double *d; // Eigenvalues
	double *e; 
} regularize_hessian2_mem;

int regularize_hessian2_calculate_mem_size(int n);
regularize_hessian_mem *regularize_hessian2_assign_mem(int n, void *raw_memory);
regularize_hessian_mem *regularize_hessian2_create_mem(int n);
int regularize_hessian2(double *H, int n, double eps, regularize_hessian_mem *mem);

#endif // REGULARIZE_HESSIAN2_H_