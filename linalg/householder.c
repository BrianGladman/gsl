#include <config.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include "gsl_linalg.h"

#include "matrix.c"

double
gsl_linalg_householder_transform (gsl_vector * v)
{
  /* replace v[0:n-1] with a householder vector (v[0:n-1]) and
     coefficient tau that annihilate v[1:n-1] */

  const size_t n = v->size ;

  double alpha, beta, tau ;

  gsl_vector x = subvector (v, 1, n-1) ; 
  
  double xnorm = gsl_blas_dnrm2 (&x);

  if (xnorm == 0) 
    {
      return 0; /* tau = 0 */
    }

  alpha = gsl_vector_get (v, 0) ;
  beta = - (alpha >= 0 ? +1 : -1) * gsl_hypot(alpha, xnorm) ;
  tau = (beta - alpha) / beta ;

  gsl_blas_dscal (1.0 / (alpha - beta), &x);
  gsl_vector_set (v, 0, beta) ;
  
  return tau;
}

int
gsl_linalg_householder_hm (double tau, const gsl_vector *v, gsl_matrix *m, gsl_vector * work)
{
  /* applies a householder transformation v,tau to matrix m */

  size_t i, j;

  if (tau == 0)
    return GSL_SUCCESS;

  /* w = M' v */

  for (i = 0; i < m->size2; i++)
    {
      double sum = gsl_matrix_get(m,0,i);  

      for (j = 1; j < m->size1; j++)  /* note, computed for v(0) = 1 above */
        {
          sum += gsl_matrix_get(m,j,i) * gsl_vector_get(v,j);
        }
      gsl_vector_set(work, i, sum);
    }

  /* M = M - v w' */

  for (j = 0; j < m->size2; j++) 
    {
      double wj = gsl_vector_get (work, j);
      double m0j = gsl_matrix_get (m, 0, j);
      gsl_matrix_set (m, 0, j, m0j - tau *  wj);
    }

  for (i = 1; i < m->size1; i++)
    {
      double vi = gsl_vector_get (v, i);

      for (j = 0; j < m->size2; j++) 
        {
          double wj = gsl_vector_get (work, j);
          double mij = gsl_matrix_get (m, i, j);
          gsl_matrix_set (m, i, j, mij - tau * vi * wj);
        }
    }

  return GSL_SUCCESS;
}

int
gsl_linalg_householder_hv (double tau, const gsl_vector *v, gsl_vector *w)
{
  /* applies a householder transformation v to vector w */
  size_t i;
  double d = 0;
 
  if (tau == 0)
    return GSL_SUCCESS ;

  /* d = v'w */

  d = gsl_vector_get(w,0);

  for (i = 1 ; i < v->size ; i++)
    {
      d += gsl_vector_get(v,i) * gsl_vector_get(w,i);
    }

  /* w = w - tau (v) (v'w) */
  
  {
    double w0 = gsl_vector_get (w,0);
    gsl_vector_set (w, 0, w0 - tau * d);
  }

  for (i = 1; i < v->size ; i++)
    {
      double wi = gsl_vector_get (w,i);
      double vi = gsl_vector_get (v,i);
      gsl_vector_set (w, i, wi - tau * vi * d);
    }

  return GSL_SUCCESS;
}

