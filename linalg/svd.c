/* Author:  G. Jungman */

#include <config.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "gsl_linalg.h"

#define REAL double

int
gsl_linalg_decomp_SV(gsl_matrix * A,
                         gsl_matrix * Q,
                         gsl_vector * S,
                         double tolerance)
{
  if(Q->size1 != A->size2 || Q->size1 != Q->size2 || S->size != A->size2) {
    return GSL_EBADLEN;
  }
  else {
    const int nrow = A->size1;
    const int ncol = A->size2;
    int i,j,k;

    /* Initialize the rotation counter and the sweep counter. */
    int count = 1;
    int sweep = 0;
    int sweepmax = ncol; 

    /* Always do at least 12 sweeps. */
    sweepmax = GSL_MAX(sweepmax, 12);

    /* Set Q to the identity matrix. */
    gsl_matrix_set_identity (Q);

    /* Orthogonalize A by plane rotations. */
    while(count > 0 && sweep <= sweepmax){

      /* Initialize rotation counter. */
      count = ncol*(ncol-1)/2;        

      for(j=0; j<ncol-1; j++){
        for(k=j+1; k<ncol; k++){
          double p = 0.0;
          double q = 0.0;
          double r = 0.0;
          double cosine, sine;
          double v;

          for(i=0;i<nrow;i++){
            /* quantities in rotation angles */
            const REAL Aij = gsl_matrix_get(A, i, j);
            const REAL Aik = gsl_matrix_get(A, i, k);
            p += Aij * Aik;
            q += Aij * Aij;
            r += Aik * Aik;
          }

          /* NOTE: this could be handled better by scaling
	   * the calculation of the inner products above.
	   * But I'm too lazy. This will have to do. [GJ]
	   */
          if(! (q*r < GSL_DBL_MAX)) {
	    /* overflow occured or will occur */
	    return GSL_EOVRFLW;
	  }
	  if(! (q*r > GSL_DBL_MIN)) {
	    /* underflow occured or will occur */
	    return GSL_EUNDRFLW;
	  }

          if(q*r == 0.0) {
            /* column elements of A are vanishingly small */
            count--;
            continue;
          }

          if((double)(p*p)/(double)(q*r) < tolerance) {
            /* columns j,k orthogonal
	     * note that p*p/(q*r) is automatically <= 1.0
	     */
            count--;
            continue;
          }

          /* calculate rotation angles */
          if(q < r) {
            cosine = 0.0;
            sine   = 1.0;
          }
          else {
            q -= r; 
            v  = sqrt(4.0*p*p + q*q);
            cosine = sqrt((v+q)/(2.0*v));
            sine   = p/(v*cosine);
          }

          /* apply rotation to A */
          for(i=0; i<nrow; i++) {
	    const REAL Aik = gsl_matrix_get(A, i, k);
            const REAL Aij = gsl_matrix_get(A, i, j);
            gsl_matrix_set(A, i, j, Aij*cosine + Aik*sine);
            gsl_matrix_set(A, i, k,  -Aij*sine + Aik*cosine);
          }

	  /* apply rotation to Q */
          for(i=0; i<ncol; i++){
            const REAL Qij = gsl_matrix_get(Q, i, j);
            const REAL Qik = gsl_matrix_get(Q, i, k);
            gsl_matrix_set(Q, i, j, Qij*cosine + Qik*sine);
            gsl_matrix_set(Q, i, k,  -Qij*sine + Qik*cosine);
          }
        }
      }

      /* Sweep completed. */
      sweep++;
    }

    /* 
     * Orthogonalization complete. Compute singular values.
     */

    if(count > 0) {
      /* reached sweep limit */
    }

    for(j=0; j<ncol; j++) {
      double q = 0.0;

      /* Calculate singular values. */
      for(i=0; i<nrow; i++){
        const REAL Aij = gsl_matrix_get(A, i, j);
        q += Aij * Aij;
      }
      gsl_vector_set(S, j, sqrt(q));

      /* Normalize vectors. */
      for(i=0; i<nrow; i++){
        const REAL Aij = gsl_matrix_get(A, i, j);
        const REAL Sj  = gsl_vector_get(S, j);
        gsl_matrix_set(A, i, j, Aij / Sj);
      }
    }

    return GSL_SUCCESS;
  }
}
