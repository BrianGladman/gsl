#include <stdio.h>
#include <unistd.h>
#include <gsl_errno.h>

#include "source.h"

int
FUNCTION(gsl_vector,read)(const int fildes, TYPE(gsl_vector) * v)
{
  const size_t n = v->size ;
  size_t bytes = read(fildes, v->data, n * sizeof(BASE)) ;

  if (bytes != n * sizeof(BASE)) 
    {
      GSL_ERROR ("read failed", GSL_EFAILED) ;
    } ;

  return 0;
} ;

int 
FUNCTION(gsl_vector,write)(const int fildes, const TYPE(gsl_vector) * v)
{
  size_t bytes = write(fildes, v->data, v->size * sizeof(BASE)) ;

  if (bytes != v->size * sizeof(BASE)) 
    {
      GSL_ERROR ("write failed", GSL_EFAILED) ;
    } ;
  
  return 0;
} ;

int
FUNCTION(gsl_vector,fprintf)(FILE * stream, const TYPE(gsl_vector) * v,
			     const char * format)
{
  size_t i ;
  const size_t n = v->size ;

  for (i = 0 ; i < n ; i++) 
    {
      int status = fprintf(stream, format, v->data[i]) ;
      
      if (status < 0) 
	{
	  GSL_ERROR ("fprintf failed", GSL_EFAILED) ;
	}

      putc('\n', stream) ;
    } ;

  return 0;
} ;

int
FUNCTION(gsl_vector,fscanf)(FILE * stream, TYPE(gsl_vector) * v)
{
  const size_t n = v->size ;
  size_t i ;

  for (i = 0 ; i < n ; i++) 
    {
      int status = fscanf(stream, IN_FORMAT, v->data + i) ;
      if (status != 1) 
	{
	  GSL_ERROR ("fscanf failed", GSL_EFAILED) ;
	}
    } ;

  return 0;
} ;
