int
FUNCTION (gsl_vector, swap) (TYPE (gsl_vector) * v, const size_t i, const size_t j)
{
  const size_t size = v->size ;
  
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
      size_t k ;
      size_t s = MULTIPLICITY * v->stride ;

      for (k = 0; k < MULTIPLICITY; k++)
        {
          ATOMIC tmp = v->data[j*s + k];
          v->data[j*s+k] = v->data[i*s + k];
          v->data[i*s+k] = tmp;
        }
    }
  
  return GSL_SUCCESS;
}

