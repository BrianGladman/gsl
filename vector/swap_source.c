int
FUNCTION (gsl_vector, swap) (TYPE (gsl_vector) * v, const size_t i, const size_t j)
{
  ATOMIC * data = v->data ;
  const size_t size = v->size ;
  const size_t stride = v->stride ;

  if (i >= size)
    {
      GSL_ERROR("first index is out of range", GSL_EINVAL);
    }

  if (j >= size)
    {
      GSL_ERROR("second index is out of range", GSL_EINVAL);
    }

  if (i != j)
    {
      const size_t s = MULTIPLICITY * stride ;
      size_t k ;

      for (k = 0; k < MULTIPLICITY; k++)
        {
          ATOMIC tmp = data[j*s + k];
          data[j*s+k] = data[i*s + k];
          data[i*s+k] = tmp;
        }
    }
  
  return GSL_SUCCESS;
}


int
FUNCTION (gsl_vector, reverse) (TYPE (gsl_vector) * v)
{
  ATOMIC * data = v->data ;
  const size_t size = v->size ;
  const size_t stride = v->stride ;

  const size_t s = MULTIPLICITY * stride ;
  
  size_t i ;

  for (i = 0 ; i < (size / 2) ; i++)
    {
      size_t j = size - i - 1 ;
      size_t k;

      for (k = 0; k < MULTIPLICITY; k++)
        {
          ATOMIC tmp = data[j*s + k];
          data[j*s+k] = data[i*s + k];
          data[i*s+k] = tmp;
        }
    }
  
  return GSL_SUCCESS;
}

