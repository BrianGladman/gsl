int
FUNCTION (gsl_vector, swap) (TYPE (gsl_vector) * v, TYPE (gsl_vector) * w)
{
  ATOMIC * d1 = v->data ;
  ATOMIC * d2 = w->data ;
  const size_t size = v->size ;
  const size_t s1 = MULTIPLICITY * v->stride ;
  const size_t s2 = MULTIPLICITY * v->stride ;
  size_t i, k ;

  if (v->size != w->size)
    {
      GSL_ERROR("vector lengths must be equal", GSL_EINVAL);
    }

  for (i = 0; i < size; i++)
    {
      for (k = 0; k < MULTIPLICITY; k++)
        {
          ATOMIC tmp = d1[i*s1 + k];
          d1[i*s1+k] = d2[i*s2 + k];
          d2[i*s2+k] = tmp;
        }
    }
  
  return GSL_SUCCESS;
}

int
FUNCTION (gsl_vector, swap_elements) (TYPE (gsl_vector) * v, const size_t i, const size_t j)
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

