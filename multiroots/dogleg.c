
int
newton_direction (const gsl_matrix * r, const gsl_vector * qtf, gsl_vector * p)
{
  int status = gsl_la_Rsolve_impl (r, qtf, p);

  for (i = 0; i < n; i++)
    {
      double pi = gsl_vector_get (p, i);
      gsl_vector_set (p, i, -pi);
    }

  return status;
}

int
gradient_direction (const gsl_matrix * r, const gsl_vector * qtf, 
                    const gsl_vector * diag, gsl_vector * g)
{
  const size_t M = r->size1 ;
  const size_t N = r->size2 ;

  for (j = 0; j < M; j++)
    {
      double sum = 0;
      double diag;

      for (i = 0; i < N; i++)
        {
          sum += gsl_matrix_get (r, i, j) * gsl_vector_get (qtf, i);
        }

      dj = gsl_vector_get(diag,j);
      gsl_vector_set (g, j, -sum/dj);
    }

  return GSL_SUCCESS;
}

int
dogleg (const gsl_matrix * r, const gsl_matrix * qtf, 
        const gsl_vector * diag, double delta, gsl_vector *p)
{
  newton_direction (r, qtf, p);

  qnorm = enorm (diag, p);

  if (qnorm <= delta)
    return GSL_SUCCESS;

  
  gradient_direction (r, qtf, diag, p

}
