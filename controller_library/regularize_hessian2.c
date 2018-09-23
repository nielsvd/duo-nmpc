#include "regularize_hessian2.h"
#include <string.h>
#include <math.h>

static double hypot2(double x, double y) { return sqrt(x * x + y * y); }

/* Symmetric Householder reduction to tridiagonal form. */
static void tred2(int n, double *V, double *d, double *e)
{
    /* This is derived from the Algol procedures tred2 by
    Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
    Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
    Fortran subroutine in EISPACK. */

    int i, j, k;
    double f, g, h, hh;
    for (j = 0; j < n; j++)
    {
        d[j] = V[(n - 1) * n + j];
    }

    /* Householder reduction to tridiagonal form. */

    for (i = n - 1; i > 0; i--)
    {
        /* Scale to avoid under/overflow. */

        double scale = 0.0;
        double h = 0.0;
        for (k = 0; k < i; k++)
        {
            scale = scale + fabs(d[k]);
        }
        if (scale == 0.0)
        {
            e[i] = d[i - 1];
            for (j = 0; j < i; j++)
            {
                d[j] = V[(i - 1) * n + j];
                V[i * n + j] = 0.0;
                V[j * n + i] = 0.0;
            }
        }
        else
        {
            /* Generate Householder vector. */

            for (k = 0; k < i; k++)
            {
                d[k] /= scale;
                h += d[k] * d[k];
            }
            f = d[i - 1];
            g = sqrt(h);
            if (f > 0)
            {
                g = -g;
            }
            e[i] = scale * g;
            h = h - f * g;
            d[i - 1] = f - g;
            for (j = 0; j < i; j++)
            {
                e[j] = 0.0;
            }

            /* Apply similarity transformation to remaining columns. */

            for (j = 0; j < i; j++)
            {
                f = d[j];
                V[j * n + i] = f;
                g = e[j] + V[j * n + j] * f;
                for (k = j + 1; k <= i - 1; k++)
                {
                    g += V[k * n + j] * d[k];
                    e[k] += V[k * n + j] * f;
                }
                e[j] = g;
            }
            f = 0.0;
            for (j = 0; j < i; j++)
            {
                e[j] /= h;
                f += e[j] * d[j];
            }
            hh = f / (h + h);
            for (j = 0; j < i; j++)
            {
                e[j] -= hh * d[j];
            }
            for (j = 0; j < i; j++)
            {
                f = d[j];
                g = e[j];
                for (k = j; k <= i - 1; k++)
                {
                    V[k * n + j] -= (f * e[k] + g * d[k]);
                }
                d[j] = V[(i - 1) * n + j];
                V[i * n + j] = 0.0;
            }
        }
        d[i] = h;
    }

    /* Accumulate transformations. */

    for (i = 0; i < n - 1; i++)
    {
        V[(n - 1) * n + i] = V[i * n + i];
        V[i * n + i] = 1.0;
        h = d[i + 1];
        if (h != 0.0)
        {
            for (k = 0; k <= i; k++)
            {
                d[k] = V[k * n + i + 1] / h;
            }
            for (j = 0; j <= i; j++)
            {
                g = 0.0;
                for (k = 0; k <= i; k++)
                {
                    g += V[k * n + i + 1] * V[k * n + j];
                }
                for (k = 0; k <= i; k++)
                {
                    V[k * n + j] -= g * d[k];
                }
            }
        }
        for (k = 0; k <= i; k++)
        {
            V[k * n + i + 1] = 0.0;
        }
    }
    for (j = 0; j < n; j++)
    {
        d[j] = V[(n - 1) * n + j];
        V[(n - 1) * n + j] = 0.0;
    }
    V[(n - 1) * n + n - 1] = 1.0;
    e[0] = 0.0;
}

/* Symmetric tridiagonal QL algorithm. */
static void tql2(int n, double *V, double *d, double *e)
{
    /*  This is derived from the Algol procedures tql2, by
    Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
    Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
    Fortran subroutine in EISPACK. */

    int i, m, l, k;
    double g, p, r, dl1, h, f, tst1, eps;
    double c, c2, c3, el1, s, s2;

    for (i = 1; i < n; i++)
    {
        e[i - 1] = e[i];
    }
    e[n - 1] = 0.0;

    f = 0.0;
    tst1 = 0.0;
    eps = pow(2.0, -52.0);
    for (l = 0; l < n; l++)
    {
        /* Find small subdiagonal element */

        tst1 = fmax(tst1, fabs(d[l]) + fabs(e[l]));
        m = l;
        while (m < n)
        {
            if (fabs(e[m]) <= eps * tst1)
            {
                break;
            }
            m++;
        }

        /* If m == l, d[l] is an eigenvalue,
        otherwise, iterate. */

        if (m > l)
        {
            int iter = 0;
            do
            {
                iter = iter + 1;
                /* Compute implicit shift */

                g = d[l];
                p = (d[l + 1] - g) / (2.0 * e[l]);
                r = hypot2(p, 1.0);
                if (p < 0)
                {
                    r = -r;
                }
                d[l] = e[l] / (p + r);
                d[l + 1] = e[l] * (p + r);
                dl1 = d[l + 1];
                h = g - d[l];
                for (i = l + 2; i < n; i++)
                {
                    d[i] -= h;
                }
                f = f + h;

                /* Implicit QL transformation. */

                p = d[m];
                c = 1.0;
                c2 = c;
                c3 = c;
                el1 = e[l + 1];
                s = 0.0;
                s2 = 0.0;
                for (i = m - 1; i >= l; i--)
                {
                    c3 = c2;
                    c2 = c;
                    s2 = s;
                    g = c * e[i];
                    h = c * p;
                    r = hypot2(p, e[i]);
                    e[i + 1] = s * r;
                    s = e[i] / r;
                    c = p / r;
                    p = c * d[i] - s * g;
                    d[i + 1] = h + s * (c * g + s * d[i]);

                    /* Accumulate transformation. */

                    for (k = 0; k < n; k++)
                    {
                        h = V[k * n + i + 1];
                        V[k * n + i + 1] = s * V[k * n + i] + c * h;
                        V[k * n + i] = c * V[k * n + i] - s * h;
                    }
                }
                p = -s * s2 * c3 * el1 * e[l] / dl1;
                e[l] = s * p;
                d[l] = c * p;

                /* Check for convergence. */

            } while (fabs(e[l]) > eps * tst1 && iter < 20); /* (Check iteration count here.) */
        }
        d[l] = d[l] + f;
        e[l] = 0.0;
    }
}

static void eigen_decomposition(double *H, int n, regularize_hessian2_mem *mem)
{
    double *V = mem->V;
    double *d = mem->d;
    double *e = mem->e; //(double *) calloc(n, sizeof(double));
    memcpy(V,H,n*n*sizeof(double));
    tred2(n, V, d, e);
    tql2(n, V, d, e);
}

static void reconstruct_A(double *H, int n, regularize_hessian2_mem *mem)
{
    double *V = mem->V;
    double *d = mem->d;
    
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j <= i; j++)
        {
            H[i * n + j] = 0.0;
            for (int k = 0; k < n; k++)
                H[i * n + j] += V[i * n + k] * d[k] * V[j * n + k];
            H[j * n + i] = H[i * n + j];
        }
    }
}

int regularize_hessian2_calculate_mem_size(int n) {
    // V: n x n
    // d: n
    // e: n
    return (sizeof(regularize_hessian2_mem) + 2*n*sizeof(double) + n*n*sizeof(double));
}

regularize_hessian2_mem *regularize_hessian2_assign_mem(int n, void *raw_memory) {
    char *c_ptr = (char *) raw_memory;

    regularize_hessian2_mem *out = (regularize_hessian2_mem *) c_ptr;
    c_ptr += sizeof(regularize_hessian2_mem);

    out->Ǜ = (double *) c_ptr;
    c_ptr += n*n*sizeof(double);

    out->d = (double *) c_ptr;
    c_ptr += n*sizeof(double);

    out->e = (double *) c_ptr;
    c_ptr += n*sizeof(double);

    assert(c_ptr == raw_memory + regularize_hessian2_calculate_mem_size(n));

    return out;
}

regularize_hessian2_mem *regularize_hessian2_create_mem(int n) {
    int size = regularize_hessian2_calculate_mem_size(n);
    void *mem = malloc(size);
    regularize_hessian2_mem *out = regularize_hessian2_assign_mem(n, mem);
    return out;
}

int regularize_hessian2(double *H, int n, double eps, regularize_hessian2_mem *mem)
{
    // set e to zero
    memset(e,0,n*sizeof(double));

    eigen_decomposition(H, n, mem);

    for (int i = 0; i < n; i++)
    {
        if (d[i] < eps)
            d[i] = eps;
    }

    reconstruct_A(H, n, mem);
}