/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include "gsl_eigen.h"


/* The eigen_sort_impl below is not very good, but it is
 * simple and self-contained. We can implement an
 * eigen_qsort_impl later by creating an array of
 * struct { val, index } and sorting those using
 * standard library qsort().
 */


int
gsl_eigen_sort_impl(gsl_vector * eval,
                    gsl_matrix * evec,
                    gsl_eigen_sort_t sort_type)
{
  if(eval == 0 || evec == 0) {
    return GSL_EFAULT;
  }
  else if(evec->size1 != evec->size2 || eval->size != evec->size1) {
    return GSL_EBADLEN;
  }
  else {
    int N = eval->size;
    int i;

    for(i=0; i<N-1; i++) {
      int j;
      int k = i;
      double tmp = gsl_vector_get(eval, k);

      /* search for something to swap */
      for(j=i+1; j<N; j++) {
        const double eval_j = gsl_vector_get(eval,j);
        const int test = (sort_type == GSL_EIGEN_SORT_VALUE ? eval_j <= tmp : fabs(eval_j) <= fabs(tmp));
        if(test) {
          k = j;
          tmp = gsl_vector_get(eval, k);
        }
      }

      if(k != i) {
        /* swap eigenvalues */
	gsl_vector_swap(eval, i, k);

        /* swap eigenvectors */
	gsl_matrix_swap_cols(evec, i, k);
      }
    }

    return GSL_SUCCESS;
  }
}
