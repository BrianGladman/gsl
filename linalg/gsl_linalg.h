/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_LINALG_H
#define GSL_LINALG_H

#include <gsl_vector.h>
#include <gsl_matrix.h>

typedef enum {
  GSL_LA_MOD_NONE=0,
  GSL_LA_MOD_TRANSPOSE=1,
  GSL_LA_MOD_CONJUGATE=2
}
gsl_la_matrix_mod_t;


typedef enum {
  GSL_LA_EIGEN_SORT_VALUE,
  GSL_LA_EIGEN_SORT_ABSVALUE
}
gsl_la_eigen_sort_t;


/* Simple implementation of matrix multiply.
 * Calculates C = A.B
 *
 * exceptions: GSL_EFAULT, GSL_EBADLEN
 */
int
gsl_la_matmult_impl(const gsl_matrix * A, const gsl_matrix * B, gsl_matrix * C);


/* Simple implementation of matrix multiply.
 * Allows transposition of either matrix, so it
 * can compute Trans(A).B or A.Trans(B) or Trans(A).Trans(B)
 *
 * exceptions: GSL_EFAULT, GSL_EBADLEN
 */
int
gsl_la_matmult_mod_impl(const gsl_matrix * A, gsl_la_matrix_mod_t modA,
                        const gsl_matrix * B, gsl_la_matrix_mod_t modB,
                        gsl_matrix * C);


/* Singular Value Decomposition
 *
 * exceptions: 
 */
int
gsl_la_decomp_SV_impl(gsl_matrix * A,
                      gsl_matrix * Q,
                      gsl_vector * S,
                      double tolerance);

/* LU Decomposition, Crout's method
 *
 * exceptions: 
 */
int
gsl_la_decomp_LU_impl(gsl_matrix * matrix,
                      gsl_vector_int * permutation,
		      int * signum);


/* Linear Solve Using LU Decomposition
 *
 * exceptions: 
 */
int
gsl_la_solve_LU_impl(const gsl_matrix     * lu_matrix,
                     const gsl_vector_int * permutation,
                     const gsl_vector     * rhs,
		     gsl_vector           * solution);

/* Linear Solve Using Householder Transformations
 *
 * exceptions: 
 */
int
gsl_la_solve_HH_impl(gsl_matrix * matrix,
                     gsl_vector * vec);


/* Eigensolve by Jacobi Method
 *
 * The data in the matrix input is destroyed.
 *
 * exceptions: 
 */
int
gsl_la_eigen_jacobi_impl(gsl_matrix * matrix,
                         gsl_vector * eval,
                         gsl_matrix * evec,
                         unsigned int max_rot, 
                         unsigned int * nrot);


/* Sort eigensystem results based on eigenvalues.
 * Sorts in order of increasing value or increasing
 * absolute value.
 *
 * exceptions: GSL_EFAULT, GSL_EBADLEN
 */
int
gsl_la_eigen_sort_impl(gsl_vector * eval,
                       gsl_matrix * evec,
                       gsl_la_eigen_sort_t sort_type);


/* Invert by Jacobi Method
 *
 * exceptions: 
 */
int
gsl_la_invert_jacobi_impl(const gsl_matrix * matrix,
                          gsl_matrix * ainv,
                          unsigned int max_rot);





#endif  /* !GSL_LINALG_H */
