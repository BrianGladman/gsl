#include <stdlib.h>

#include <gsl_errno.h>

typedef struct
{
  size_t n;
  double * data ;
} gsl_vector ;

double
gsl_vector_get(gsl_vector v, size_t i)
{
#ifdef GSL_RANGE_CHECK
  if (i < 0 || i >= v.n) 
    {
      abort() ;
    }
#endif
  return v.data[i] ;
}


double *
gsl_vector_set(gsl_vector v, size_t i)
{
#ifdef GSL_RANGE_CHECK
  if (i < 0 || i >= v.n) 
    {
      abort() ;
    }
#endif
  return v.data + i ;
}


void
gsl_vector_set_directly(gsl_vector v, size_t i, double x)
{
#ifdef GSL_RANGE_CHECK
  if (i < 0 || i >= v.n) 
    {
      abort() ;
    }
#endif
  v.data[i] = x ;
}

int
gsl_vector_alloc (gsl_vector * v, size_t n)
{
  if (n == 0)
    {
      GSL_ERROR ("vector length n must be positive integer", GSL_EDOM);
    }
  
  v->data = malloc(n * sizeof(double)) ;

  if (v->data == 0) 
    {
      GSL_ERROR ("failed to allocate space for vector", GSL_ENOMEM);
    }
  
  v->n = n ;

  return 0 ;
}

int
gsl_vector_calloc (gsl_vector * v, size_t n)
{
  size_t i ;

  int status = gsl_vector_alloc (v, n) ;
  
  if (status) 
    return status ;

  for (i = 0 ; i < n; i++)  /* initialize to zero */
    {
      v->data[i] = 0.0 ;
    }

  return 0 ;
}


int
gsl_vector_free (gsl_vector * v)
{
  free(v->data) ;
  v->n = 0 ;
  return 0 ;
}




