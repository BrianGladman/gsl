void
generate_householder (gsl_vector * v, double * tau, size_t i)
{
  /* replace v[0:n-1] with a householder vector (tau, v[1:n-1]) that
     annihilates v[1:n-1] */

  double x0 = gsl_vector_get (v, 0) ;

  gsl_vector v1 = subvector (v, 1, n-1) ; 
  
  double s = gsl_blas_dnrm2 (&v1);
  
  double mu;

  if (s == 0) 
    {
      *tau = 0 ;
      return ;
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

  *tau  =  (2 * v0) / (s * s + v0 * v0) ;

  gsl_blas_dscal (1.0 / v0, &c);
  
  return ;
}

void
apply_householder (gsl_matrix *m, double tau, gsl_vector * work, size_t i)
{
  /* applies a householder transformation in column i to matrix m */
 
  double mii = gsl_matrix_get (m, i, i) ;
  gsl_matrix_set (m, i, i, 1) ;

  gsl_matrix a = submatrix (m, i+1, i+1, m, n);

  gsl_vector v = col (m, i, i, m) ;

  /* w = M' v */

  gsl_blas_dgemv (CblasTrans, 1.0, &a, &v, 0.0, work) ;

  /* M = M - w v' */

  gsl_blas_dger (-tau, work, &v, &a);

}
