BASE 
FUNCTION(gsl_matrix,max) (const TYPE(gsl_matrix) * m)
{
  /* finds the largest element of a matrix */

  const size_t M = m->size1 ;
  const size_t N = m->size2 ;

  BASE max = FUNCTION(gsl_matrix,get)(m,0,0);
  size_t i,j;

  for (i = 0; i < M; i++)
    {
      for (j = 0; j < N; j++)
        {
          BASE x = FUNCTION(gsl_matrix,get)(m, i, j);
          if (x > max)
            max = x;
        }
    }

  return max;
}

BASE 
FUNCTION(gsl_matrix,min) (const TYPE(gsl_matrix) * m)
{
  /* finds the smallest element of a matrix */

  const size_t M = m->size1 ;
  const size_t N = m->size2 ;

  BASE min = FUNCTION(gsl_matrix,get)(m,0,0);
  size_t i,j;

  for (i = 0; i < M; i++)
    {
      for (j = 0; j < N; j++)
        {
          BASE x = FUNCTION(gsl_matrix,get)(m, i, j);
          if (x < min)
            min = x;
        }
    }

  return min;
}


void
FUNCTION(gsl_matrix,minmax) (const TYPE(gsl_matrix) * m,
                             BASE * min_out, size_t * imin, size_t * jmin, 
                             BASE * max_out, size_t * imax, size_t * jmax)
{
  /* finds the smallest and largest elements of a matrix */

  const size_t M = m->size1 ;
  const size_t N = m->size2 ;

  BASE max = FUNCTION(gsl_matrix,get)(m,0,0);
  BASE min = FUNCTION(gsl_matrix,get)(m,0,0);
  size_t i,j;

  for (i = 0; i < M; i++)
    {
      for (j = 0; j < N; j++)
        {
          BASE x = FUNCTION(gsl_matrix,get)(m, i, j);
          if (x < min)
            {
              min = x;
              *imin = i;
              *jmin = j;
            }
          if (x > max)
            {
              max = x;
              *imax = i;
              *jmax = j;
            }
        }
    }

  *min_out = min;
  *max_out = max;
}
