int
FUNCTION (gsl_matrix, memcpy) (TYPE (gsl_matrix) * dest,
                               const TYPE (gsl_matrix) * src)
{
  const size_t src_size1 = src->size1;
  const size_t src_size2 = src->size2;
  const size_t dest_size1 = dest->size1;
  const size_t dest_size2 = dest->size2;

  if (src_size1 != dest_size1 || src_size2 != dest_size2)
    {
      GSL_ERROR ("matrix sizes are different", GSL_EBADLEN);
    }

  {
    const size_t src_dim2 = src->dim2 ;
    const size_t dest_dim2 = dest->dim2 ;
    size_t i, j;

    for (i = 0; i < src_size1 ; i++)
      {
        for (j = 0; j < MULTIPLICITY * src_size2; j++)
          {
            dest->data[dest_dim2 * i + j] = src->data[src_dim2 * i + j];
          }
      }
  }

  return GSL_SUCCESS;
}

