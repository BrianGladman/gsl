void
generate_householder (gsl_matrix * m, gsl_vector * tau, size_t i)
{
  /* replace [tau(i), m(i+1:m,i)] with a householder vector that
     annihilates m(i:m,i) */
  
  sigma = column_norm (m, i+1, n, i) ;
  
  if (sigma == 0) 
    {
      gsl_vector_set (tau, i, 0);
      return ;
    }

  x0 = gsl_matrix_get (m, i, i) ;

  mu = sqrt(x0 * x0 + sigma) ;

  if (x0 <= 0)
    {
      v0 = x0 - mu ;
    }
  else
    {
      v0 = -sigma / (x0 + mu) ;
    }

  gsl_vector_set (tau, i, 2 * v0 / (sigma + v0 * v0)) ;

  column_scale (m, i+1, n, i, 1.0/v0) ;
}

void
apply_householder (gsl_matrix *m, double tau, gsl_vector * work, size_t i)
{
  /* applies a householder transformation in column i to matrix m */
 
  /* w = M' v */

  double mii = gsl_matrix_get (m, i, i) ;

  gsl_matrix_set (m, i, i, 1) ;

  gsl_blas_dgemv (CblasTrans, 1.0, submatrix(m, i+1, i+1, m, n),
                  col(m, i, i, m), 0.0, work) ;

  gsl_blas_dger (-tau, work, col(m, i, i, m), submatrix(m, i+1, i+1, m, n));

}
