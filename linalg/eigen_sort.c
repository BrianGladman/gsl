/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <stdlib.h>
#include <gsl_math.h>
#include "gsl_linalg.h"


/* The eigen_sort_impl below is not very good, but it is
 * simple and self-contained. We can implement an
 * eigen_qsort_impl later by creating an array of
 * struct { val, index } and sorting those using
 * standard library qsort().
 */


int
gsl_la_eigen_sort_impl(gsl_vector * eval,
                       gsl_matrix * evec,
                       gsl_la_eigen_sort_t sort_type)
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
    double tmp;
    gsl_vector * tmp_vec_1 = gsl_vector_alloc(N);
    gsl_vector * tmp_vec_2 = gsl_vector_alloc(N);

    for(i=0; i<N-1; i++) {
      int j;
      int k = i;
      tmp = eval->data[k];

      /* search for something to swap */
      for(j=i+1; j<N; j++) {
        int test = (sort_type == GSL_LA_EIGEN_SORT_VALUE ? eval->data[j] <= tmp : fabs(eval->data[j]) <= fabs(tmp));
        if(test) {
          k = j;
          tmp = eval->data[k];
        }
      }

      if(k != i) {
        /* swap eigenvalues */
        eval->data[k] = eval->data[i];
        eval->data[i] = tmp;

        /* swap eigenvectors */ /* matrix should probably export row/col swap ops */
        gsl_matrix_copy_col(evec, i, tmp_vec_1);
	gsl_matrix_copy_col(evec, k, tmp_vec_2);
	gsl_matrix_set_col(evec, i, tmp_vec_2);
	gsl_matrix_set_col(evec, k, tmp_vec_1);
      }
    }

    gsl_vector_free(tmp_vec_1);
    gsl_vector_free(tmp_vec_2);
    return GSL_SUCCESS;
  }
}
