/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_EIGEN_H
#define GSL_EIGEN_H

#include <gsl_vector.h>
#include <gsl_matrix.h>

typedef enum {
  GSL_EIGEN_SORT_VALUE,
  GSL_EIGEN_SORT_ABSVALUE
}
gsl_eigen_sort_t;


/* Eigensolve by Jacobi Method
 *
 * The data in the matrix input is destroyed.
 *
 * exceptions: 
 */
int
gsl_eigen_jacobi_impl(gsl_matrix * matrix,
                      gsl_vector * eval,
                      gsl_matrix * evec,
                      unsigned int max_rot, 
                      unsigned int * nrot);


/* Invert by Jacobi Method
 *
 * exceptions: 
 */
int
gsl_eigen_invert_jacobi_impl(const gsl_matrix * matrix,
                             gsl_matrix * ainv,
                             unsigned int max_rot);


/* Sort eigensystem results based on eigenvalues.
 * Sorts in order of increasing value or increasing
 * absolute value.
 *
 * exceptions: GSL_EFAULT, GSL_EBADLEN
 */
int
gsl_eigen_sort_impl(gsl_vector * eval,
                    gsl_matrix * evec,
                    gsl_eigen_sort_t sort_type);


#endif  /* !GSL_EIGEN_H */
