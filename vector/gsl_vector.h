
typedef struct
{
  size_t n;
  size_t n_allocated;
  double * data ;
} gsl_vector ;

double
gsl_vector_get(gsl_vector V, size_t i)
{
#ifdef GSL_RANGE_CHECK
  if (i < 0 || i >= n) 
    {
      abort() ;
    }
#endif
  return V.data[i] ;
}


double *
gsl_vector_set(gsl_vector V, size_t i)
{
#ifdef GSL_RANGE_CHECK
  if (i < 0 || i >= n) 
    {
      abort() ;
    }
#endif
  return V.data + i ;
}


int 
gsl_vector_alloc (
