/* eigen/eigen_sort.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

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
	gsl_vector_swap_elements(eval, i, k);

        /* swap eigenvectors */
	gsl_matrix_swap_columns(evec, i, k);
      }
    }

    return GSL_SUCCESS;
  }
}
