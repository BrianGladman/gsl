#include <gsl_errno.h>

int
FUNCTION(gsl_matrix,get_row)(const TYPE(gsl_matrix) * m,
			     const size_t i,
			     TYPE (gsl_vector) * v)
{
  const size_t column_length = m->size1 ;
  const size_t row_length = m->size2 ;

  if (i >= column_length)
    {
      GSL_ERROR("row index is out of range", GSL_EINVAL) ;
    }

  if (v->size != row_length) 
    {
      GSL_ERROR("matrix row size and vector length are not equal", 
		GSL_EBADLEN) ;
    }

  
  v->data = m->data + MULTIPLICITY * i * row_length;
  v->stride = 1;

  return 0 ;
}

int
FUNCTION(gsl_matrix,get_col)(const TYPE(gsl_matrix) * m,
			     const size_t j,
			     TYPE (gsl_vector) * v)
{
  const size_t column_length = m->size1 ;
  const size_t row_length = m->size2 ;

  if (j >= row_length)
    {
      GSL_ERROR("column index is out of range", GSL_EINVAL) ;
    }

  if (v->size != column_length) 
    {
      GSL_ERROR("matrix column size and vector length are not equal", 
		GSL_EBADLEN) ;
    }

  
  v->data = m->data + MULTIPLICITY * j;
  v->stride = row_length;

  return 0 ;
}


int
FUNCTION(gsl_matrix,set_row)(TYPE(gsl_matrix) * m,
			     const size_t i,
			     const TYPE (gsl_vector) * v)
{
  const size_t column_length = m->size1 ;
  const size_t row_length = m->size2 ;

  if (i >= column_length)
    {
      GSL_ERROR("row index is out of range", GSL_EINVAL) ;
    }

  if (v->size != row_length) 
    {
      GSL_ERROR("matrix row size and vector length are not equal", 
		GSL_EBADLEN) ;
    }

  {
    const BASE * v_data = (BASE *) v->data ;
    BASE * row_data = (BASE *) m->data + MULTIPLICITY * i * row_length ;
    size_t j ;

    for (j = 0 ; j < row_length ; j++) 
      {
	row_data[j] = v_data[j] ;
      }
  }

  return 0 ;
}

int
FUNCTION(gsl_matrix,set_col)(TYPE(gsl_matrix) * m,
			     const size_t j,
			     const TYPE (gsl_vector) * v)
{
  const size_t column_length = m->size1 ;
  const size_t row_length = m->size2 ;

  if (j >= row_length)
    {
      GSL_ERROR("column index is out of range", GSL_EINVAL) ;
    }

  if (v->size != column_length) 
    {
      GSL_ERROR("matrix column size and vector length are not equal", 
		GSL_EBADLEN) ;
    }

  {
    const BASE * v_data = (BASE *)v->data ;
    BASE * column_data = ((BASE *)m->data) + j ;
    size_t i ;

    for (i = 0 ; i < column_length ; i++) 
      {
	column_data[i*row_length] = v_data[i] ;
      }
  }

  return 0 ;
}


