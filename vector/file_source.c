#include <stdio.h>
#include <unistd.h>
#include <gsl_errno.h>

#include "source.h"

int
FUNCTION(gsl_vector,fread)(FILE * stream, TYPE(gsl_vector) * v)
{
  int status = FUNCTION(gsl_block,fread) (stream, v->data, v->size) ;
  return status;
}

int 
FUNCTION(gsl_vector,fwrite)(FILE * stream, const TYPE(gsl_vector) * v)
{
  int status = FUNCTION(gsl_block,fwrite) (stream, v->data, v->size) ;
  return status;
}

int
FUNCTION(gsl_vector,fprintf)(FILE * stream, const TYPE(gsl_vector) * v,
			     const char * format)
{
  int status = FUNCTION(gsl_block,fprintf) (stream, v->data, v->size, format) ;
  return status;
}

int
FUNCTION(gsl_vector,fscanf)(FILE * stream, TYPE(gsl_vector) * v)
{
  int status = FUNCTION(gsl_block,fscanf) (stream, v->data, v->size) ;
  return status;
}

int
FUNCTION(gsl_block,fread)(FILE * stream, BASE * data, const size_t n)
{
  size_t items = fread (data,  sizeof(BASE), n, stream) ;

  if (items != n) 
    {
      GSL_ERROR ("read failed", GSL_EFAILED) ;
    } ;

  return 0;
}

int 
FUNCTION(gsl_block,fwrite)(FILE * stream, const BASE * data, const size_t n)
{
  size_t items = fwrite(data, sizeof(BASE), n, stream) ;

  if (items != n) 
    {
      GSL_ERROR ("write failed", GSL_EFAILED) ;
    } ;
  
  return 0;
}


int
FUNCTION(gsl_block,fprintf)(FILE * stream, const BASE * data, const size_t n,
			     const char * format)
{
  size_t i ;

  for (i = 0 ; i < n ; i++) 
    {
      int status = fprintf(stream, format, data[i]) ;
      
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
      int status = fscanf(stream, IN_FORMAT, data + i) ;

      if (status != 1) 
	{
	  GSL_ERROR ("fscanf failed", GSL_EFAILED) ;
	}
    }

  return 0;
} 


