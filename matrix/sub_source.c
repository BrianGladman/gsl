#include <gsl_errno.h>

int
FUNCTION(gsl_matrix,get_row)(const TYPE(gsl_matrix) * m,
			     const size_t i,
			     TYPE (gsl_vector) * v)
{
  const size_t column_range = m->size1 ;
  const size_t row_length = m->size2 ;

  if (i >= column_range)
    {
      GSL_ERROR("row index is out of range", GSL_EINVAL) ;
    }

  if (v->size != row_length) 
    {
      GSL_ERROR("matrix row size and vector length are not equal", 
		GSL_EBADLEN) ;
    }

  {
    const BASE * v_data = v->data ;
    const BASE * row_data = m->data + i * row_length ;
    size_t j ;

    for (j = 0 ; j < row_length ; j++) 
      {
	v_data[j] = row_data[j] ;
      }
  }

  return 0 ;
}

int
FUNCTION(gsl_matrix,set_row)(TYPE(gsl_matrix) * m,
			     const size_t i,
			     const TYPE (gsl_vector) * v)
{
  const size_t column_range = m->size1 ;
  const size_t row_length = m->size2 ;

  if (i >= column_range)
    {
      GSL_ERROR("row index is out of range", GSL_EINVAL) ;
    }

  if (v->size != row_length) 
    {
      GSL_ERROR("matrix row size and vector length are not equal", 
		GSL_EBADLEN) ;
    }

  {
    const BASE * v_data = v->data ;
    const BASE * row_data = m->data + i * row_length ;
    size_t j ;

    for (j = 0 ; j < row_length ; j++) 
      {
	row_data[j] = v_data[j] ;
      }
  }

  return 0 ;
}


