static double column_norm (gsl_matrix * a, size_t i0, size_t i1, size_t j);

static double column_norm (gsl_matrix * a, size_t i0, size_t i1, size_t j)
{
  size_t i;
  REAL sum = 0;

  for (i = i0; i < i1; i++)
    {
      REAL aij = gsl_matrix_get(a, i, j);

      sum += aij * aij;   /* FIXME: use the proper scaled algorithm here */
    }
  return sqrt(sum);
}
