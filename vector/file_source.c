#include <stdio.h>
#include <unistd.h>
#include <gsl_errno.h>


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
FUNCTION(gsl_block,fread)(FILE * stream, ATOMIC * data, const size_t n)
{
  size_t items = fread (data,  sizeof(BASE), n, stream) ;

  if (items != n) 
    {
      GSL_ERROR ("fread failed", GSL_EFAILED) ;
    } ;

  return 0;
}

int 
FUNCTION(gsl_block,fwrite)(FILE * stream, const ATOMIC * data, const size_t n)
{
  size_t items = fwrite(data, sizeof(BASE), n, stream) ;

  if (items != n) 
    {
      GSL_ERROR ("fwrite failed", GSL_EFAILED) ;
    } ;
  
  return 0;
}


