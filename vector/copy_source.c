int
FUNCTION (gsl_vector, copy) (TYPE (gsl_vector) * dest,
                             const TYPE (gsl_vector) * src)
{
  const size_t src_size = src->size;
  const size_t dest_size = dest->size;

  if (src_size != dest_size)
    {
      GSL_ERROR ("vector lengths are not equal", GSL_EBADLEN);
    }

  {
    const size_t src_stride = src->stride ;
    const size_t dest_stride = dest->stride ;
    size_t j;

    for (j = 0; j < src_size; j++)
      {
        size_t k;

        for (k = 0; k < MULTIPLICITY; k++) 
          {
            dest->data[MULTIPLICITY * dest_stride * j + k] 
              = src->data[MULTIPLICITY * src_stride*j + k];
          }
      }
  }

  return GSL_SUCCESS;
}

