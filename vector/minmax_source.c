BASE 
FUNCTION(gsl_vector,max) (const TYPE(gsl_vector) * v)
{
  /* finds the largest element of a vector */

  const size_t N = v->size ;
  const size_t stride = v->stride ;

  BASE max = v->data[0 * stride];
  size_t i;

  for (i = 0; i < N; i++)
    {
      BASE x = v->data[i*stride];
      if (x > max)
        max = x;
    }

  return max;
}

BASE 
FUNCTION(gsl_vector,min) (const TYPE(gsl_vector) * v)
{
  /* finds the smallest element of a vector */

  const size_t N = v->size ;
  const size_t stride = v->stride ;

  BASE min = v->data[0 * stride];
  size_t i;

  for (i = 0; i < N; i++)
    {
      BASE x = v->data[i*stride];
      if (x < min)
        min = x;
    }

  return min;
}

void
FUNCTION(gsl_vector,minmax) (const TYPE(gsl_vector) * v,
                             BASE * min_out, 
                             BASE * max_out)
{
  /* finds the smallest and largest elements of a vector */

  const size_t N = v->size ;
  const size_t stride = v->stride ;

  BASE max = v->data[0 * stride];
  BASE min = v->data[0 * stride];

  size_t i;

  for (i = 0; i < N; i++)
    {
      BASE x = v->data[i*stride];
      if (x < min)
        {
          min = x;
        }
      if (x > max)
        {
          max = x;
        }
    }

  *min_out = min;
  *max_out = max;
}


size_t 
FUNCTION(gsl_vector,max_index) (const TYPE(gsl_vector) * v)
{
  /* finds the largest element of a vector */

  const size_t N = v->size ;
  const size_t stride = v->stride ;

  BASE max = v->data[0 * stride];
  size_t imax = 0;
  size_t i;

  for (i = 0; i < N; i++)
    {
      BASE x = v->data[i*stride];
      if (x > max)
        {
          max = x;
          imax = i;
        }
    }

  return imax;
}

size_t 
FUNCTION(gsl_vector,min_index) (const TYPE(gsl_vector) * v)
{
  /* finds the smallest element of a vector */

  const size_t N = v->size ;
  const size_t stride = v->stride ;

  BASE min = v->data[0 * stride];
  size_t imin = 0;
  size_t i;

  for (i = 0; i < N; i++)
    {
      BASE x = v->data[i*stride];
      if (x < min)
        {
          min = x;
          imin = i;
        }
    }

  return imin;
}


void
FUNCTION(gsl_vector,minmax_index) (const TYPE(gsl_vector) * v,
                                   size_t * imin_out, 
                                   size_t * imax_out)
{
  /* finds the smallest and largest elements of a vector */

  const size_t N = v->size ;
  const size_t stride = v->stride ;

  size_t imin = 0, imax = 0;
  BASE max = v->data[0 * stride];
  BASE min = v->data[0 * stride];

  size_t i;

  for (i = 0; i < N; i++)
    {
      BASE x = v->data[i*stride];
      if (x < min)
        {
          min = x;
          imin = i;
        }
      if (x > max)
        {
          max = x;
          imax = i;
        }
    }

  *imin_out = imin;
  *imax_out = imax;
}


