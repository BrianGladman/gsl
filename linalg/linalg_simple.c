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
#include <gsl_math.h>
#include <gsl_vector.h>
#include <gsl_matrix.h>
#include "gsl_linalg.h"

#define REAL double


int
gsl_la_decomp_SV_impl(gsl_matrix * A,
                      gsl_matrix * Q,
                      gsl_vector * S,
                      double tolerance)
{
  if(A == 0 || Q == 0 || S == 0) {
    return GSL_EFAULT;
  }
  else if(Q->size1 != A->size2 || Q->size1 != Q->size2 || S->size != A->size2) {
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
    for(i=0; i<ncol; i++){
      for(j=0; j<ncol; j++){
        gsl_matrix_set(Q, i, j, 0.0);
      }
      gsl_matrix_set(Q, i, i, 1.0);
    }

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

          if(q*r < GSL_DBL_EPSILON) {
            /* column elements of A small */
            count--;
            continue;
          }

          if(p*p/(q*r) < tolerance) {
            /* columns j,k orthogonal */
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
	    REAL Aik = gsl_matrix_get(A, i, k);
            REAL Aij = gsl_matrix_get(A, i, j);
            gsl_matrix_set(A, i, j, Aij*cosine + Aik*sine);
            gsl_matrix_set(A, i, k,  -Aij*sine + Aik*cosine);
          }
          for(i=0; i<ncol; i++){     /* apply rotation to Q */
            REAL Qij = gsl_matrix_get(Q, i, j);
            REAL Qik = gsl_matrix_get(Q, i, k);
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
        REAL Aij = gsl_matrix_get(A, i, j);
        q += Aij * Aij;
      }
      gsl_vector_set(S, j, sqrt(q));

      /* Normalize vectors. */
      for(i=0; i<nrow; i++){
        REAL Aij = gsl_matrix_get(A, i, j);
        REAL Sj  = gsl_vector_get(S, j);
        gsl_matrix_set(A, i, j, Aij / Sj);
      }
    }

    return GSL_SUCCESS;
  }
}


int
gsl_la_decomp_LU_impl(gsl_matrix * matrix,
                      gsl_vector_int * permutation,
		      int * signum)
{
  if(matrix == 0 || permutation == 0) {
    return GSL_EFAULT;
  }
  else if(matrix->size1 != matrix->size2) {
    return GSL_ENOTSQR;
  }
  else if(permutation->size != matrix->size1) {
    return GSL_EBADLEN;
  }
  else if(matrix->size1 == 0) {
    return GSL_SUCCESS; /* FIXME: what to do with this idiot case ? */
  }
  else {
    int i, j;
    const int N  = matrix->size1;
    REAL * scale = (REAL *) malloc(N * sizeof(REAL));

    if(scale == 0) {
      return GSL_ENOMEM;
    }

    /* Prepare permutation and scaling information.
     */
    *signum = 1;
    for(i=0; i<N; i++) {
      REAL max_row_element = 0.0;
      gsl_vector_int_set(permutation, i, i);

      for(j=0; j<N; j++) {
        REAL aij = fabs(gsl_matrix_get(matrix, i, j));
        max_row_element = GSL_MAX(max_row_element, aij);
      }

      if(max_row_element == 0.0) {
        /* Trap exact singularity.
	 */
        *signum = 0;
        free(scale);
        return GSL_ESING;
      }

      scale[i] = 1.0/max_row_element;
    }

    /* Perform decomposition.
     */
    for(i=0; i<N; i++) {

      int perm_index_i = gsl_vector_int_get(permutation, i);
      REAL scale_i = scale[perm_index_i];
      REAL elem_ij = gsl_matrix_get(matrix, perm_index_i, j);
      REAL pivot   = fabs(elem_ij) * scale_i;
      int  j_pivot = i;

      for(j=i+1; j<N; j++) {
        int perm_index_j = gsl_vector_int_get(permutation, j);
	REAL scale_j = scale[perm_index_j];
        REAL tmp = fabs(gsl_matrix_get(matrix, perm_index_j, i)) * scale_j;
        if(pivot < tmp) {
          pivot   = tmp;
          j_pivot = j;
        }
      }

      if(pivot < GSL_DBL_EPSILON) {
        /* Trap approximate singularity.
	 */
        *signum = 0;
        free(scale);
        return GSL_ESING;
      }

      /* Execute the pivot (if it is non-null).
       */
      if(j_pivot != i) {
        int perm_tmp1 = gsl_vector_int_get(permutation, j_pivot);
        int perm_tmp2 = gsl_vector_int_get(permutation, i);
	gsl_vector_int_set(permutation, i, perm_tmp1);
	gsl_vector_int_set(permutation, j_pivot, perm_tmp2);
        *signum = - *signum;
      }

      /* Execute Gauss step.
       */
      for(j=i+1; j<N; j++) {
        int perm_index_j = gsl_vector_int_get(permutation, j);
	REAL elem_ji = gsl_matrix_get(matrix, perm_index_j, i);
	REAL elem_ii = gsl_matrix_get(matrix, perm_index_i, i);
        if(elem_ji != 0.0) {
	  int m;
	  REAL tmp;
	  gsl_matrix_set(matrix, perm_index_j, i, elem_ji / elem_ii);
	  tmp = gsl_matrix_get(matrix, perm_index_j, i);
	  for(m=i+1; m<N; m++) {
	    REAL elem_im = gsl_matrix_get(matrix, perm_index_i, m);
	    REAL elem_jm = gsl_matrix_get(matrix, perm_index_j, m);
	    gsl_matrix_set(matrix, perm_index_j, m, elem_jm - tmp * elem_im);
          }
        }
      }
    }

    free(scale);
    return GSL_SUCCESS;
  }
}


int
gsl_la_solve_LU_impl(const gsl_matrix     * lu_matrix,
                     const gsl_vector_int * permutation,
                     const gsl_vector     * rhs,
		     gsl_vector           * solution)
{
  if(lu_matrix == 0 || permutation == 0 || rhs == 0 || solution == 0) {
    return GSL_EFAULT;
  }
  else if(lu_matrix->size1 != lu_matrix->size2) {
    return GSL_ENOTSQR;
  }
  else if(   lu_matrix->size1 != permutation->size
          || lu_matrix->size1 != rhs->size
	  || lu_matrix->size1 != solution->size)
         {
    return GSL_EBADLEN;
  }
  else if(lu_matrix->size1 == 0) {
    return GSL_SUCCESS; /* FIXME: dumb case */
  }
  else {
    int j, k;
    const int N = lu_matrix->size1;

    /* Apply permutation to RHS
     * and perform update.
     */
    for(k=0; k<N; k++) {
      int perm_index_k = gsl_vector_int_get(permutation, k);
      gsl_vector_set(solution, k, gsl_vector_get(rhs, perm_index_k));
      for(j=0; j<k; j++) {
        REAL sol_k  = gsl_vector_get(solution, k);
	REAL sol_j  = gsl_vector_get(solution, j);
	REAL lum_kj = gsl_matrix_get(lu_matrix, k, j);
        gsl_vector_set(solution, k, sol_k - lum_kj * sol_j);
      }
    }

    /* Perform back-substitution.
     */
    for(k=N-1; k>=0; k--) {
      REAL lum_kk;
      REAL sum = 0.0;
      for(j=k+1; j<N; j++) {
        REAL lum_kj = gsl_matrix_get(lu_matrix, k, j);
	REAL sol_j  = gsl_vector_get(solution, j);
        sum += lum_kj * sol_j;
      }

      lum_kk = gsl_matrix_get(lu_matrix, k, k);
      if(lum_kk == 0.0) {
        return GSL_EINVAL;
      }
      else {
        REAL sol_k = gsl_vector_get(solution, k);
	gsl_vector_set(solution, k, (sol_k - sum) / lum_kk);
      }
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



#ifdef HAVE_INLINE
inline
#endif
static void
jac_rotate(gsl_matrix * a,
           int i, int j, int k, int l,
           double * g, double * h,
           double s, double tau)
{
  *g = gsl_matrix_get(a, i, j);
  *h = gsl_matrix_get(a, k, l);
  gsl_matrix_set(a, i, j, (*g) - s*((*h) + (*g)*tau));
  gsl_matrix_set(a, k, l, (*h) + s*((*g) - (*h)*tau));
}


int
gsl_la_eigen_jacobi_impl(gsl_matrix * a,
                         gsl_vector * eval,
                         gsl_matrix * evec,
                         unsigned int max_rot, 
                         unsigned int * nrot)
{
  if(a == 0 || eval == 0 || evec == 0) {
    return GSL_EFAULT;
  }
  else if(a->size1 != a->size2) {
    return GSL_ENOTSQR;
  }
  else if(a->size1 != evec->size1 || a->size1 != evec->size2) {
    return GSL_EBADLEN;
  }
  else if(a->size1 != eval->size) {
    return GSL_EBADLEN;
  }
  else {
    const int n = a->size1;
    unsigned int i;
    int j, iq, ip;
    double t, s;

    REAL * b = (REAL *) malloc(n * sizeof(REAL));
    REAL * z = (REAL *) malloc(n * sizeof(REAL));
    if(b == 0 || z == 0) {
      if(b != 0) free(b);
      if(z != 0) free(z);
      return GSL_ENOMEM;
    }

    /* Set eigenvectors to coordinate basis. */
    for(ip=0; ip<n; ip++) {
      for(iq=0; iq<n; iq++) {
        gsl_matrix_set(evec, ip, iq, 0.0);
      }
      gsl_matrix_set(evec, ip, ip, 1.0);
    }

    /* Initialize eigenvalues and workspace. */
    for(ip=0; ip<n; ip++) {
      REAL a_ipip = gsl_matrix_get(a, ip, ip);
      z[ip] = 0.0;
      b[ip] = a_ipip;
      gsl_vector_set(eval, ip, a_ipip);
    }

    *nrot = 0;

    for(i=1; i<=max_rot; i++) {
      REAL thresh;
      REAL tau;
      REAL g, h, c;
      REAL sm = 0.0;
      for(ip=0; ip<n-1; ip++) {
        for(iq=ip+1; iq<n; iq++) {
          sm += fabs(gsl_matrix_get(a, ip, iq));
        }
      }
      if(sm == 0.0) {
        free(z);
        free(b);
        return GSL_SUCCESS;
      }

      if(i < 4)
        thresh = 0.2*sm/(n*n);
      else
        thresh = 0.0;

      for(ip=0; ip<n-1; ip++) {
        for(iq=ip+1; iq<n; iq++) {
          const REAL d_ip = gsl_vector_get(eval, ip);
          const REAL d_iq = gsl_vector_get(eval, iq);
	  const REAL a_ipiq = gsl_matrix_get(a, ip, iq);
          g = 100.0 * fabs(a_ipiq);
          if(   i > 4
             && fabs(d_ip)+g == fabs(d_ip)
             && fabs(d_iq)+g == fabs(d_iq)
	     ) {
            gsl_matrix_set(a, ip, iq, 0.0);
          }
          else if(fabs(a_ipiq) > thresh) {
            h = d_iq - d_ip;
            if(fabs(h) + g == fabs(h)) {
              t = a_ipiq/h;
            }
            else {
              REAL theta = 0.5*h/a_ipiq;
              t = 1.0/(fabs(theta) + sqrt(1.0 + theta*theta));
              if(theta < 0.0) t = -t;
            }

            c   = 1.0/sqrt(1.0+t*t);
            s   = t*c;
            tau = s/(1.0+c);
            h   = t * a_ipiq;
            z[ip] -= h;
            z[iq] += h;
	    gsl_vector_set(eval, ip, d_ip - h);
	    gsl_vector_set(eval, iq, d_iq + h);
	    gsl_matrix_set(a, ip, iq, 0.0);

            for(j=0; j<=ip-1; j++){
              jac_rotate(a, j, ip, j, iq, &g, &h, s, tau);
            }
            for(j=ip+1; j<=iq-1; j++){
              jac_rotate(a, ip, j, j, iq, &g, &h, s, tau);
            }
            for(j=iq+1; j<n; j++){
              jac_rotate(a, ip, j, iq, j, &g, &h, s, tau);
            }
            for (j=0; j<n; j++){
              jac_rotate(evec, j, ip, j, iq, &g, &h, s, tau);
            }
            ++(*nrot);
          }
        }
      }
      for (ip=0; ip<n; ip++) {
        b[ip] += z[ip];
        z[ip]  = 0.0;
	gsl_vector_set(eval, ip, b[ip]);
      }

      /* continue iteration */
    }

    return GSL_EMAXITER;
  }
}


int
gsl_la_invert_jacobi_impl(const gsl_matrix * a,
                          gsl_matrix * ainv,
                          unsigned int max_rot)
{
  if(a == 0 || ainv == 0) {
    return GSL_EFAULT;
  }
  else if(a->size1 != a->size2 || ainv->size1 != ainv->size2) {
    return GSL_ENOTSQR;
  }
  else if(a->size1 != ainv->size2) {
    return GSL_EBADLEN;
  }
  else {
    const int n = a->size1;
    unsigned int nrot;
    int i,j,k,l;

    /* This is annoying because I do not want
     * the error handling in these functions.
     * But there are no "impl"-like versions
     * of these allocators... sigh.
     */
    gsl_vector * eval = gsl_vector_alloc(n);
    gsl_matrix * evec = gsl_matrix_alloc(n, n);
    gsl_matrix * inv_diag = gsl_matrix_alloc(n, n);

    if(eval == 0 || evec == 0 || inv_diag == 0) {
      if(eval != 0) gsl_vector_free(eval);
      if(evec != 0) gsl_matrix_free(evec);
      if(inv_diag != 0) gsl_matrix_free(inv_diag);
      return GSL_ENOMEM;
    }

    memcpy(ainv->data, a->data, n*n*sizeof(REAL));

    gsl_la_eigen_jacobi_impl(ainv, eval, evec, max_rot, &nrot);

    for(i=0; i<n; i++) {
      if(fabs(gsl_vector_get(eval, i)) < GSL_SQRT_DBL_EPSILON) {
        /* apparent singularity */
        gsl_vector_free(eval);
        gsl_matrix_free(evec);
        gsl_matrix_free(inv_diag);
        return GSL_ESING;
      }
    }

    /* Invert the diagonalized matrix. */
    for(i=0; i<n; i++) {
      for(j=0; j<n; j++) {
        gsl_matrix_set(inv_diag, i, j, 0.0);
      }
      gsl_matrix_set(inv_diag, i, i, 1.0/gsl_vector_get(eval, i));
    }

    for(i=0; i<n; i++) {
      for(j=0; j<n; j++) {
        gsl_matrix_set(ainv, i, j, 0.0);
        for(k=0; k<n; k++) {
          for(l=0; l<n; l++) {
	    REAL ainv_ij = gsl_matrix_get(ainv, i, j);
	    REAL evec_il = gsl_matrix_get(evec, i, l);
	    REAL evec_jk = gsl_matrix_get(evec, j, k);
	    REAL inv_diag_lk = gsl_matrix_get(inv_diag, l, k);
	    REAL delta = evec_il * inv_diag_lk * evec_jk;
	    gsl_matrix_set(ainv, i, j, ainv_ij + delta);
	  }
	}
      }
    }

    gsl_vector_free(eval);
    gsl_matrix_free(evec);
    gsl_matrix_free(inv_diag);
    return GSL_SUCCESS;
  }
}
