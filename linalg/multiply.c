/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <stdlib.h>
#include <gsl_math.h>
#include "gsl_linalg.h"



int
gsl_la_matmult_impl(const gsl_matrix * A, const gsl_matrix * B, gsl_matrix * C)
{
  if(A == 0 || B == 0 || C == 0) {
    return GSL_EFAULT;
  }
  else if(A->size2 != B->size1 || A->size1 != C->size1 || B->size2 != C->size2) {
    return GSL_EBADLEN;
  }
  else {
    double a, b;
    double temp;
    size_t i, j, k;

    for(i=0; i<C->size1; i++) {
      for(j=0; j<C->size2; j++) {
        a = A->data[i*A->size2 + 0];
	b = B->data[0 + j];
        temp = a * b;
        for(k=1; k<A->size2; k++) {
	  a = A->data[i*A->size2 + k];
	  b = B->data[k*B->size2 + j];
          temp += a * b;
	}
	C->data[i*C->size2 + j] = temp;
      }
    }

    return GSL_SUCCESS;
  }
}
