/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef __GSL_EIGEN_H__
#define __GSL_EIGEN_H__

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

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


__END_DECLS

#endif /* __GSL_EIGEN_H__ */
