#include <config.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_permutation.h>

size_t
gsl_permutation_size (const gsl_permutation * p)
{
  return p->size ;
}

size_t *
gsl_permutation_data (const gsl_permutation * p)
{
  return p->data ;
}

size_t
gsl_permutation_get (const gsl_permutation * p, const size_t i)
{
  if (gsl_check_range)
    {
      if (i >= p->size)		/* size_t is unsigned, can't be negative */
	{
	  GSL_ERROR_RETURN ("index out of range", GSL_EINVAL, 0);
	}
    }

  return p->data[i];
}

int
gsl_permutation_swap (gsl_permutation * p, const size_t i, const size_t j)
{
  const size_t size = p->size ;
  
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
      size_t tmp = p->data[i];
      p->data[i] = p->data[j];
      p->data[j] = tmp;
    }
  
  return GSL_SUCCESS;
}


int
gsl_permutation_valid (gsl_permutation * p)
{
  const size_t size = p->size ;

  size_t i, j ;

  for (i = 0; i < size; i++) 
    {
      if (p->data[i] >= size)
        {
          GSL_ERROR("permutation index outside range", GSL_FAILURE) ;
        }

      for (j = 0; j < i; j++)
        {
          if (p->data[i] == p->data[j])
            {
              GSL_ERROR("duplicate permutation index", GSL_FAILURE) ;
            }
        }
    }
  
  return GSL_SUCCESS;
}



int
gsl_permutation_reverse (gsl_permutation * p)
{
  const size_t size = p->size ;

  size_t i ;
  
  for (i = 0; i < size; i++) 
    {
      size_t tmp ;
      size_t j = size - i - 1;

      if (i == j)
        continue ;
      
      tmp = p->data[i] ;
      p->data[i] = p->data[j] ;
      p->data[j] = tmp ;
    }
}
