/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <stdlib.h>
#include <gsl_math.h>
#include "gsl_linalg.h"


#define SWAP_SIZE_T(a, b)  do { size_t swap_tmp = a; a = b; b = swap_tmp; } while(0)


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


int
gsl_la_matmult_mod_impl(const gsl_matrix * A, gsl_la_matrix_mod_t modA,
                        const gsl_matrix * B, gsl_la_matrix_mod_t modB,
                        gsl_matrix * C)
{
  if(A == 0 || B == 0 || C == 0) {
    return GSL_EFAULT;
  }
  else if(modA == GSL_LA_MOD_NONE && modB == GSL_LA_MOD_NONE) {
    return gsl_la_matmult_impl(A, B, C);
  }
  else {
    size_t dim1_A = A->size1;
    size_t dim2_A = A->size2;
    size_t dim1_B = B->size1;
    size_t dim2_B = B->size2;
    size_t dim1_C = C->size1;
    size_t dim2_C = C->size2;

    if(modA & GSL_LA_MOD_TRANSPOSE) SWAP_SIZE_T(dim1_A, dim2_A);
    if(modB & GSL_LA_MOD_TRANSPOSE) SWAP_SIZE_T(dim1_B, dim2_B);

    if(dim2_A != dim1_B || dim1_A != dim1_C || dim2_B != dim2_C) {
      return GSL_EBADLEN;
    }
    else {
      double a, b;
      double temp;
      size_t i, j, k;
      size_t a1, a2, b1, b2;

      for(i=0; i<dim1_C; i++) {
        for(j=0; j<dim2_C; j++) {
          a1 = i;
          a2 = 0;
          b1 = 0;
          b2 = j;
          if(modA & GSL_LA_MOD_TRANSPOSE) SWAP_SIZE_T(a1, a2);
	  if(modB & GSL_LA_MOD_TRANSPOSE) SWAP_SIZE_T(b1, b2);

          a = A->data[a1*dim2_A + a2];
          b = B->data[b1*dim2_B + b2];
          temp = a * b;

          for(k=1; k<dim2_A; k++) {
	    a1 = i;
            a2 = k;
            b1 = k;
            b2 = j;
            if(modA & GSL_LA_MOD_TRANSPOSE) SWAP_SIZE_T(a1, a2);
	    if(modB & GSL_LA_MOD_TRANSPOSE) SWAP_SIZE_T(b1, b2);
            a = A->data[a1*dim2_A + a2];
            b = B->data[b1*dim2_B + b2];
            temp += a * b;
          }

          C->data[i*dim2_C + j] = temp;
	}
      }

      return GSL_SUCCESS;
    }
  }
}
