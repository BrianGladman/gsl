double
generate_householder (gsl_vector * v)
{
  /* replace v[0:n-1] with a householder vector (tau, v[1:n-1]) that
     annihilates v[1:n-1] */

  const size_t n = v->size ;

  double x0 = gsl_vector_get (v, 0) ;

  gsl_vector v1 = subvector (v, 1, n-1) ; 
  
  double s = gsl_blas_dnrm2 (&v1);
  
  double mu, v0;

  if (s == 0) 
    {
      return 0; /* tau = 0 */
    }

  mu = gsl_hypot(x0, s) ;

  if (x0 <= 0)
    {
      v0 = x0 - mu ;
    }
  else
    {
      v0 = -(s*s) / (x0 + mu) ;
    }

  gsl_blas_dscal (1.0 / v0, &v1);

  gsl_vector_set (v, 0, mu) ;
  
  return (2 * v0) / (s * s + v0 * v0);  /* this is tau */
}

void
apply_householder (gsl_vector *v, gsl_matrix *m, double tau, gsl_vector * work)
{
  /* applies a householder transformation in column i to matrix m */

  double v0;
 
  if (tau == 0)
    return ;

  v0 = gsl_vector_get (v, 0);

  gsl_vector_set (v, 0, 1.0) ;

  /* w = M' v */

  gsl_blas_dgemv (CblasTrans, 1.0, m, v, 0.0, work) ;

  /* M = M - v w' */

  gsl_blas_dger (-tau, v, work, m);

  gsl_vector_set (v, 0, v0);
}

void
apply_householder_v (gsl_vector *v, gsl_vector *w, double tau)
{
  /* applies a householder transformation to vector w */

  double v0,d;
 
  if (tau == 0)
    return ;

  v0 = gsl_vector_get (v, 0);

  gsl_vector_set (v, 0, 1.0) ;

  /* d = v'w */

  gsl_blas_ddot (v, w, &d);

  /* w = w - tau (v) (v'w) */

  gsl_blas_daxpy (-tau * d, v, w);

  gsl_vector_set (v, 0, v0);
}
