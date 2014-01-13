#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sparse.h>

/* exact solution */
double u_exact(const double x) { return sin(M_PI * x); }

int
main()
{
  const size_t N = 100;                       /* number of grid points */
  const size_t n = N - 2;                     /* subtract 2 to exclude boundaries */
  const double h = 1.0 / (N - 1.0);           /* grid spacing */
  gsl_spmatrix *T = gsl_spmatrix_alloc(n ,n); /* triplet format */
  gsl_spmatrix *C;                            /* compressed format */
  gsl_vector *f = gsl_vector_alloc(n);        /* right hand side vector */
  gsl_vector *u = gsl_vector_alloc(n);        /* solution vector */
  size_t i;

  /* construct the sparse matrix for the finite difference equation */

  /* loop over interior grid points */
  for (i = 0; i < n; ++i)
    {
      /* u_{i+1} term, ignore at the boundary */
      if (i + 1 < n)
        gsl_spmatrix_set(T, i, i + 1, 1.0);

      gsl_spmatrix_set(T, i, i, -2.0);

      /* u_{i-1} term, ignore at the boundary */
      if (i > 0)
        gsl_spmatrix_set(T, i, i - 1, 1.0);
    }

  /* scale by h^2 */
  gsl_spmatrix_scale(T, 1.0 / (h * h));

  /* construct right hand side vector */
  for (i = 0; i < n; ++i)
    {
      double xi = (i + 1) * h;
      double fi = -M_PI * M_PI * sin(M_PI * xi);
      gsl_vector_set(f, i, fi);
    }

  /* convert to compressed column format */
  C = gsl_spmatrix_compress(T);

  /*
   * At this point, C->i, C->p, C->data contain the
   * row indices, column pointers, and matrix elements
   * of the compressed column storage format. These 3
   * arrays can be passed to external linear solvers.
   *
   * For illustration purposes, we will convert the
   * sparse matrix to a dense gsl_matrix and complete
   * the solution using a dense LU solver.
   */

  {
    gsl_matrix *A = gsl_matrix_alloc(n, n);
    gsl_permutation *p = gsl_permutation_alloc(n);
    int s;

    /* convert sparse to dense */
    gsl_spmatrix_sp2d(A, T);

    /* solve linear system A u = f */
    gsl_linalg_LU_decomp(A, p, &s);
    gsl_linalg_LU_solve(A, p, f, u);

    /* output solution */
    for (i = 0; i < n; ++i)
      {
        double xi = (i + 1) * h;
        double u_analytic = u_exact(xi);
        double u_gsl = gsl_vector_get(u, i);

        printf("%f %.12e %.12e\n", xi, u_gsl, u_analytic);
      }

    gsl_matrix_free(A);
    gsl_permutation_free(p);
  }

  gsl_spmatrix_free(T);
  gsl_spmatrix_free(C);
  gsl_vector_free(f);
  gsl_vector_free(u);

  return 0;
} /* main() */
