/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_LINALG_H
#define GSL_LINALG_H

#include <gsl_vector.h>
#include <gsl_matrix.h>


/* Singular Value Decomposition
 *
 * exceptions: 
 */
int
gsl_la_decomp_SV_impl(gsl_matrix * A,
                      gsl_matrix * Q,
                      gsl_vector * S,
                      double tolerance);

/* LU Decomposition
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
 * exceptions: 
 */
int
gsl_la_eigen_jacobi_impl(gsl_matrix * a,
                         gsl_vector * eval,
                         gsl_matrix * evec,
                         unsigned int max_rot, 
                         unsigned int * nrot);

/* Invert by Jacobi Method
 */
int
gsl_la_invert_jacobi_impl(const gsl_matrix * a,
                          gsl_matrix * ainv,
                          unsigned int max_rot);


#endif  /* !GSL_LINALG_H */
