int
FUNCTION (gsl_matrix, copy_row) (TYPE (gsl_vector) * v,
                                 const TYPE (gsl_matrix) * m,
				 const size_t i)
{
  const size_t column_length = m->size1;
  const size_t row_length = m->size2;

  if (i >= column_length)
    {
      GSL_ERROR ("row index is out of range", GSL_EINVAL);
    }

  if (v->size != row_length)
    {
      GSL_ERROR ("matrix row size and vector length are not equal",
		 GSL_EBADLEN);
    }

  {
    ATOMIC *row_data = m->data + MULTIPLICITY * i * row_length;
    const size_t stride = v->stride ;
    size_t j;

    for (j = 0; j < MULTIPLICITY * row_length; j++)
      {
	v->data[stride * j] = row_data[j];
      }
  }

  return GSL_SUCCESS;
}

int
FUNCTION (gsl_matrix, copy_col) (TYPE (gsl_vector) * v,
                                 const TYPE (gsl_matrix) * m,
				 const size_t j)
{
  const size_t column_length = m->size1;
  const size_t row_length = m->size2;

  if (j >= row_length)
    {
      GSL_ERROR ("column index is out of range", GSL_EINVAL);
    }

  if (v->size != column_length)
    {
      GSL_ERROR ("matrix column size and vector length are not equal",
		 GSL_EBADLEN);
    }


  {
    ATOMIC *column_data = m->data + MULTIPLICITY * j;
    const size_t stride = v->stride ;
    size_t i;

    for (i = 0; i < row_length; i++)
      {
	size_t k;
	for (k = 0; k < MULTIPLICITY; k++)
	  {
	    v->data[stride * MULTIPLICITY * j + k] =
	      column_data[MULTIPLICITY * i * row_length + k];
	  }
      }
  }

  return GSL_SUCCESS;
}

int
FUNCTION (gsl_matrix, set_row) (TYPE (gsl_matrix) * m,
				const size_t i,
				const TYPE (gsl_vector) * v)
{
  const size_t column_length = m->size1;
  const size_t row_length = m->size2;

  if (i >= column_length)
    {
      GSL_ERROR ("row index is out of range", GSL_EINVAL);
    }

  if (v->size != row_length)
    {
      GSL_ERROR ("matrix row size and vector length are not equal",
		 GSL_EBADLEN);
    }

  {
    const ATOMIC *v_data = v->data;
    ATOMIC *row_data = m->data + MULTIPLICITY * i * row_length;
    const size_t stride = v->stride ;
    size_t j;

    for (j = 0; j < MULTIPLICITY * row_length; j++)
      {
	row_data[j] = v_data[stride * j];
      }
  }

  return GSL_SUCCESS;
}

int
FUNCTION (gsl_matrix, set_col) (TYPE (gsl_matrix) * m,
				const size_t j,
				const TYPE (gsl_vector) * v)
{
  const size_t column_length = m->size1;
  const size_t row_length = m->size2;

  if (j >= row_length)
    {
      GSL_ERROR ("column index is out of range", GSL_EINVAL);
    }

  if (v->size != column_length)
    {
      GSL_ERROR ("matrix column size and vector length are not equal",
		 GSL_EBADLEN);
    }

  {
    const ATOMIC *v_data = v->data;
    ATOMIC *column_data = m->data + MULTIPLICITY * j;
    const size_t stride = v->stride ;
    size_t i;

    for (i = 0; i < column_length; i++)
      {
	column_data[MULTIPLICITY * i * row_length] = v_data[stride * i];
      }
  }

  return GSL_SUCCESS;
}


TYPE (gsl_vector) *
FUNCTION (gsl_vector, alloc_row_from_matrix) (TYPE(gsl_matrix) * m,
                                              const size_t i)
{
  TYPE (gsl_vector) * v;

  const size_t column_length = m->size1;

  if (i >= column_length)
    {
      GSL_ERROR_RETURN ("row index is out of range", GSL_EINVAL, 0);
    }

  v = (TYPE (gsl_vector) *) malloc (sizeof (TYPE (gsl_vector)));

  if (v == 0)
    {
      GSL_ERROR_RETURN ("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
    }

  v->data = m->data + i * m-> dim2 ;
  v->size = m->size2;
  v->stride = 1;
  v->block = 0;

  return v;
}

TYPE (gsl_vector) *
FUNCTION (gsl_vector, alloc_col_from_matrix) (TYPE(gsl_matrix) * m,
                                              const size_t j)
{
  TYPE (gsl_vector) * v;

  const size_t row_length = m->size2;

  if (j >= row_length)
    {
      GSL_ERROR_RETURN ("column index is out of range", GSL_EINVAL, 0);
    }

  v = (TYPE (gsl_vector) *) malloc (sizeof (TYPE (gsl_vector)));

  if (v == 0)
    {
      GSL_ERROR_RETURN ("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
    }

  v->data = m->data + j ;
  v->size = m->size1;
  v->stride = m->dim2;
  v->block = 0;

  return v;
}
