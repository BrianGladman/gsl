/* Author:  G. Jungman
 * RCS:     $Id$
 */
/* Simple linear algebra operations, operating directly
 * on the gsl_vector and gsl_matrix objects. These are
 * meant for "generic" and "small" systems. Anyone
 * interested in large systems will want to use more
 * sophisticated methods, presumably involving native
 * BLAS operations, specialized data representations,
 * or other optimizations.
 */
#include <config.h>
#include <stdlib.h>
#include <string.h>
#include <gsl_math.h>
#include <gsl_vector.h>
#include <gsl_matrix.h>
#include "gsl_linalg.h"

#define REAL double

#include "norm.c"

/* From SLATEC, qrfac */

int
gsl_la_decomp_QR_impl(gsl_matrix * matrix, 
                      gsl_vector * rdiag,
                      gsl_vector_int * permutation, 
                      int * signum)
{
  if(matrix == 0 || permutation == 0 || signum == 0) {
    return GSL_EFAULT;
  }
  else if(matrix->size1 != matrix->size2) {
    return GSL_ENOTSQR;
  }
  else if(permutation->size != matrix->size2) {
    return GSL_EBADLEN;
  }
  else if(rdiag->size != matrix->size2) {
    return GSL_EBADLEN;
  }
  else if(matrix->size1 == 0) {
    return GSL_SUCCESS; /* FIXME: what to do with this idiot case ? */
  }
  else {
    const int N = matrix->size1;
    const int M = matrix->size2;
    int i, j, k;
    REAL * wa = (REAL *) malloc(N * sizeof(REAL));

    if(wa == 0) {
      return GSL_ENOMEM;
    }

    *signum = 1;

    for (j = 0; j < N; j++)
      {
        REAL norm = column_norm(matrix, 0, M, j);
        gsl_vector_set(rdiag,j,norm);
        wa[j] = norm;
        gsl_vector_int_set(permutation,j,j) ;
      }

    for (j = 0; j < N; j++)
      {
        REAL ajnorm, rk = 0, rkmax = 0;

        /* Bring the column of largest norm into the pivot position */

        int kmax = j;

        for (k = j ; k < N; k++)
          {
            rk = gsl_vector_get(rdiag, k);

            if (rk > rkmax) 
              {
                kmax = k;
                rkmax = rk;
              }
          }

        if (kmax != k)
          {
            gsl_matrix_swap_cols (matrix, j, kmax);
            gsl_vector_int_swap (permutation, j, kmax);
            gsl_vector_set(rdiag, kmax, gsl_vector_get(rdiag,j));
            wa[kmax] = wa[j];
            (*signum) =  - (*signum);
          }

        /* Compute the Householder transformation to reduce the j-th
           column of the matrix to a multiple of the j-th unit vector */

        ajnorm = column_norm(matrix, j, M, j);

        if (ajnorm == 0)
          {
            gsl_vector_set(rdiag,j,0.0); 
            continue ;
          }
        
        if (gsl_matrix_get(matrix, j, j) < 0)
          ajnorm *= -1;

        for (i = j; i < M; i++)
          {
            REAL aij = gsl_matrix_get(matrix, i, j);
            gsl_matrix_set(matrix, i, j, aij/ajnorm);
          }

        gsl_matrix_set(matrix, j, j, 1.0 + gsl_matrix_get(matrix, j, j));

        /* Apply the transformation to the remaining columns and
           update the norms */

        for (k = j + 1; k < N ; k++)
          {
            REAL temp, sum = 0.0;

            for (i = j; i < M; i++)
              {
                REAL aij = gsl_matrix_get(matrix, i, j);
                REAL aik = gsl_matrix_get(matrix, i, k);
                sum += aij * aik;
              }

            temp = sum / gsl_matrix_get(matrix, j, j);
            
            for (i = j; i < M; i++)
              {
                REAL aij = gsl_matrix_get(matrix, i, j);
                REAL aik = gsl_matrix_get(matrix, i, k);
                gsl_matrix_set(matrix, i, k, aik - temp * aij);
              }
            
            rk = gsl_vector_get(rdiag,k); 

            if (rk == 0)
              continue;
            
            temp = gsl_matrix_get(matrix, j, k) / rk ;

            if (fabs(temp) >= 1)
              rk = 0.0;
            else
              rk = rk*sqrt(1-temp*temp);
            
            if (fabs(rk/wa[k]) < sqrt(20.0) * GSL_SQRT_DBL_EPSILON)
              {
                rk = column_norm(matrix, j+1, M, k);
                wa[k] = rk;
              }

            gsl_vector_set(rdiag,k, rk);
          }

        gsl_vector_set(rdiag,j,-ajnorm);
      }
    free(wa);
    return GSL_SUCCESS;
  }
}

int
gsl_la_solve_QR_impl(const gsl_matrix     * qr_matrix,
                     const gsl_vector     * rdiag,
                     const gsl_vector_int * permutation,
                     const gsl_vector     * rhs,
		     gsl_vector           * solution)
{
  if(qr_matrix == 0 || permutation == 0 || rhs == 0 || solution == 0) {
    return GSL_EFAULT;
  }
  else if(solution->data == 0 || rhs->data == 0) {
    return GSL_EFAULT;
  }
  else if(qr_matrix->size1 != qr_matrix->size2) {
    return GSL_ENOTSQR;
  }
  else if(   qr_matrix->size1 != permutation->size
          || qr_matrix->size1 != rhs->size
	  || qr_matrix->size1 != solution->size)
         {
    return GSL_EBADLEN;
  }
  else if(qr_matrix->size1 == 0) {
    return GSL_SUCCESS; /* FIXME: dumb case */
  }
  else {
    const size_t N = qr_matrix->size1;
    size_t i,j,n;
    int k;

    /* Copy rhs into solution */

    for(k=0; k<N; k++) {
      gsl_vector_set(solution, k, gsl_vector_get(rhs, k));
    }

    for (j = 0; j < N; j++)
      { 
        /* compute Q^T b */

        REAL t, sum = 0.0;

        for (i = j; i < N; i++)
          {
            REAL qij = gsl_matrix_get(qr_matrix, i, j);
            REAL bi = gsl_vector_get(solution, i);

            sum += qij * bi ;
          }

        t = -sum / gsl_matrix_get(qr_matrix, j, j);

        for (i = j; i < N; i++)
          {
            REAL qij = gsl_matrix_get(qr_matrix, i, j);
            REAL bi = gsl_vector_get(solution, i);

            gsl_vector_set(solution, i, bi + t * qij) ;
          }
      }

    /* Solve R z = Q^T b */

    for (n = N ; n > 0; n--)
      {
        size_t j, i = n - 1;

        REAL bi, di, sum = 0.0;

        for (j = n; j < N; j++)
          {
            REAL aij = gsl_matrix_get(qr_matrix, i, j);
            REAL bj = gsl_vector_get(solution, j);
            sum += aij * bj;
          }
        
        bi = gsl_vector_get(solution, i);
        di = gsl_vector_get(rdiag, i);

        gsl_vector_set(solution, i, (bi - sum)/di);
      }

    /* Apply permutation to solution */

    for(k=0; k<N; k++) {
      int perm_index_k = gsl_vector_int_get(permutation, k);
      gsl_vector_set(solution, perm_index_k, gsl_vector_get(solution, k));
    }

    return GSL_SUCCESS;
  }
}


/* [Engeln-Mullges + Uhlig, Alg. 4.42]
 */
int
gsl_la_solve_HH_impl(gsl_matrix * matrix,
                     gsl_vector * vec)
{
  if(matrix == 0 || vec == 0) {
    return GSL_EFAULT;
  }
  else if(matrix->size1 > matrix->size2) {
    /* System is underdetermined.
     */
    return GSL_EINVAL;
  }
  else if(matrix->size2 != vec->size) {
    return GSL_EBADLEN;
  }
  else if(matrix->size1 == 0 || matrix->size2 == 0) {
    return GSL_SUCCESS; /* FIXME: dumb case */
  }
  else {
    const int N = matrix->size1;
    const int M = matrix->size2;
    int i, j, k;
    REAL * d = (REAL *) malloc(N * sizeof(REAL));

    if(d == 0) {
      return GSL_ENOMEM;
    }

    /* Perform Householder transformation.
     */
    for(i=0; i<N; i++) {
      const REAL elem_ii = gsl_matrix_get(matrix, i, i);
      REAL alpha;
      REAL f;
      REAL ak;
      REAL max_norm = 0.0;
      REAL r = 0.0;

      for(k=i; k<M; k++) {
        REAL elem_ki = gsl_matrix_get(matrix, k, i);
        r += elem_ki * elem_ki;
      }

      if(r == 0.0) {
        /* Rank of matrix is
	 * less than size1.
	 */
        free(d);
        return GSL_ESING;
      }

      alpha = sqrt(r) * GSL_SIGN(elem_ii);

      ak = 1.0 / (r + alpha * elem_ii);
      gsl_matrix_set(matrix, i, i, elem_ii + alpha);

      d[i] = -alpha;

      for(k=i+1; k<N; k++) {
        REAL norm = 0.0;
	f = 0.0;
        for(j=i; j<M; j++) {
	  REAL elem_jk = gsl_matrix_get(matrix, j, k);
	  REAL elem_ji = gsl_matrix_get(matrix, j, i);
	  norm += elem_jk * elem_jk;
	  f    += elem_jk * elem_ji;
        }
	max_norm = GSL_MAX(max_norm, norm);

        f *= ak;

        for(j=i; j<M; j++) {
	  REAL elem_jk = gsl_matrix_get(matrix, j, k);
	  REAL elem_ji = gsl_matrix_get(matrix, j, i);
	  gsl_matrix_set(matrix, j, k, elem_jk - f * elem_ji);
	}
      }

      if(fabs(alpha) < 2.0 * GSL_DBL_EPSILON * sqrt(max_norm)) {
        /* Apparent singularity.
         */
        free(d);
        return GSL_ESING;
      }

      /* Perform update of RHS.
       */
      f = 0.0;
      for(j=i; j<M; j++) {
        f += gsl_vector_get(vec, j) * gsl_matrix_get(matrix, j, i);
      }
      f *= ak;
      for(j=i; j<M; j++) {
	REAL vec_j = gsl_vector_get(vec, j);
	gsl_vector_set(vec, j, vec_j - f * gsl_matrix_get(matrix, j, i));
      }

    }

    /* Perform back-substitution.
     */
    for(i=N-1; i>=0; i--) {
      REAL vec_i = gsl_vector_get(vec, i);
      REAL sum = 0.0;
      for(k=i+1; k<N; k++) {
        sum += gsl_matrix_get(matrix, i, k) * gsl_vector_get(vec, k);
      }

      gsl_vector_set(vec, i, (vec_i - sum) / d[i]);
    }

    free(d);
    return GSL_SUCCESS;
  }
}

int
gsl_la_qform_QR_impl (const gsl_matrix * qr, gsl_matrix * q)
{
  if(qr == 0 || q == 0) {
    return GSL_EFAULT;
  }
  else if(qr->size1 != qr->size2) {
    return GSL_ENOTSQR;
  }
  else if(q->size1 != qr->size2) {
    return GSL_EBADLEN;
  }
  else if(qr->size1 == 0) {
    return GSL_SUCCESS; /* FIXME: what to do with this idiot case ? */
  }
  else {
    const int N = qr->size1;
    const int M = qr->size2;
    int i, j, k, L;
    REAL * wa = (REAL *) malloc(N * sizeof(REAL));

    if(wa == 0) {
      return GSL_ENOMEM;
    }
    
    for (j = 0; j < M; j++)
      {
        for (i = 0; i < j; i++)
          gsl_matrix_set(q, i, j, 0.0);

        /* gsl_matrix_set(q, j, j, 1.0); */
        
        for (i = j ; i < M; i++)
          gsl_matrix_set(q, i, j, gsl_matrix_get(qr, i, j));
      }
    
    for (L = 0; L < M; L++)
      {
        k = (M-1) - L;

        for (i = k ; i < M; i++)
          {
            wa[i] = gsl_matrix_get(q,i,k);
            gsl_matrix_set(q,i,k,0.0);
          }

        gsl_matrix_set(q,k,k,1.0);
        
        if (wa[k] == 0)
          continue;

        for (j = k; j < M; j++)
          {
            double temp, sum = 0.0;
            for (i = k; i < M; i++)
              {
                sum += gsl_matrix_get(q,i,j) * wa[i];
              }
            temp = sum / wa[k];
            for (i = k; i < M; i++)
              {
                REAL qij = gsl_matrix_get(q, i, j);
                gsl_matrix_set(q, i, j, qij - temp*wa[i]);
              }
          }
      }
    free(wa);
    return GSL_SUCCESS;
  }
}

