#ifndef REGULARIZE_HESSIAN_H_
#define REGULARIZE_HESSIAN_H_

typedef struct regularize_hessian_mem_ {
	double *V; // Eigenvectors
	double *d; // Eigenvalues
	double *e; 
} regularize_hessian_mem;

int regularize_hessian_calculate_mem_size(int n);
regularize_hessian_mem *regularize_hessian_assign_mem(int n, void *raw_memory);
regularize_hessian_mem *regularize_hessian_create_mem(int n);
int regularize_hessian(double *H, int n, double eps, regularize_hessian_mem *mem);

#endif // REGULARIZE_HESSIAN_H_