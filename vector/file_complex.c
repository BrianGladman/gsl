#include <config.h>
#include <gsl_vector_complex.h>

#define BASE gsl_complex
#define SHORT complex
#define IN_FORMAT "%lg %lg"

#include "file_source.c"

int
FUNCTION(gsl_block,fprintf)(FILE * stream, const BASE * data, const size_t n,
			     const char * format)
{
  size_t i ;

  for (i = 0 ; i < n ; i++) 
    {
      int status = fprintf(stream, format, data[i].real) ;
      
      if (status < 0) 
	{
	  GSL_ERROR ("fprintf failed", GSL_EFAILED) ;
	}

      status = putc(' ', stream) ;

      if (status == EOF) 
	{
	  GSL_ERROR ("putc failed", GSL_EFAILED) ;
	}

      status = fprintf(stream, format, data[i].imag) ;
      
      if (status < 0) 
	{
	  GSL_ERROR ("fprintf failed", GSL_EFAILED) ;
	}

      status = putc('\n', stream) ;

      if (status == EOF) 
	{
	  GSL_ERROR ("putc failed", GSL_EFAILED) ;
	}
    }

  return 0 ;
}

int
FUNCTION(gsl_block,fscanf)(FILE * stream, BASE * data, const size_t n)
{
  size_t i ;

  for (i = 0 ; i < n ; i++) 
    {
      double x, y ;
      int status = fscanf(stream, IN_FORMAT, &x, &y) ;

      (data + i)->real = x ;
      (data + i)->imag = y ;

      if (status != 2) 
	{
	  GSL_ERROR ("fscanf failed", GSL_EFAILED) ;
	}
    }

  return 0;
} 



