/**********************************************************************/
/* The functions below are obsolete                                   */
/**********************************************************************/

int
FUNCTION (gsl_matrix, get_row) (TYPE (gsl_vector) * v,
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
FUNCTION (gsl_matrix, get_col) (TYPE (gsl_vector) * v,
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
      GSL_ERROR_VAL ("row index is out of range", GSL_EINVAL, 0);
    }

  v = (TYPE (gsl_vector) *) malloc (sizeof (TYPE (gsl_vector)));

  if (v == 0)
    {
      GSL_ERROR_VAL ("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
    }

  v->data = m->data + MULTIPLICITY * i * m->tda ;
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
      GSL_ERROR_VAL ("column index is out of range", GSL_EINVAL, 0);
    }

  v = (TYPE (gsl_vector) *) malloc (sizeof (TYPE (gsl_vector)));

  if (v == 0)
    {
      GSL_ERROR_VAL ("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
    }

  v->data = m->data + MULTIPLICITY * j ;
  v->size = m->size1;
  v->stride = m->tda;
  v->block = 0;

  return v;
}

int
FUNCTION (gsl_vector, view_row_from_matrix) (TYPE(gsl_vector) * v,
                                             TYPE(gsl_matrix) * m,
                                             const size_t i)
{
  const size_t column_length = m->size1;

  if (i >= column_length)
    {
      GSL_ERROR ("row index is out of range", GSL_EINVAL);
    }

  if (v->block != 0)
    {
      GSL_ERROR ("vector already has memory allocated to it", GSL_ENOMEM);
    }

  v->data = m->data + MULTIPLICITY * i * m->tda ;
  v->size = m->size2;
  v->stride = 1;

  return GSL_SUCCESS;
}

int
FUNCTION (gsl_vector, view_col_from_matrix) (TYPE(gsl_vector) * v,
                                             TYPE(gsl_matrix) * m,
                                             const size_t j)
{
  const size_t row_length = m->size2;

  if (j >= row_length)
    {
      GSL_ERROR_VAL ("column index is out of range", GSL_EINVAL, 0);
    }

  if (v->block != 0)
    {
      GSL_ERROR ("vector already has memory allocated to it", GSL_ENOMEM);
    }

  v->data = m->data + MULTIPLICITY * j ;
  v->size = m->size1;
  v->stride = m->tda;

  return GSL_SUCCESS;
}

int
FUNCTION(gsl_matrix, view_from_vector) (TYPE(gsl_matrix) * m,
                                        TYPE(gsl_vector) * base,
                                        size_t offset, 
                                        size_t n1, size_t n2, size_t d2)
{
  if (n1 == 0)
    {
      GSL_ERROR ("matrix dimension n1 must be positive integer", GSL_EDOM);
    }
  else if (n2 == 0)
    {
      GSL_ERROR ("matrix dimension n2 must be positive integer", GSL_EDOM);
    }
  else if (d2 < n2)
    {
      GSL_ERROR ("matrix dimension d2 must be greater than n2",	GSL_EDOM);
    }
  else if (base->size < offset + n1 * d2)
    {
      GSL_ERROR ("matrix size exceeds available vector size", GSL_EDOM);
    }
  else if (base->stride != 1) 
    {
      GSL_ERROR ("matrix view requires a vector with unit stride", GSL_EDOM);
    }

  if (m->block != 0)
    {
      GSL_ERROR ("matrix already has memory allocated to it", GSL_ENOMEM);
    }

  m->data = base->data + MULTIPLICITY * offset;
  m->size1 = n1;
  m->size2 = n2;
  m->tda = d2;
  
  return GSL_SUCCESS;
}
