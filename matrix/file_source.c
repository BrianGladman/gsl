#include <stdio.h>
#include <unistd.h>
#include <gsl_errno.h>

int
FUNCTION(gsl_matrix,fread)(FILE * stream, TYPE(gsl_matrix) * m)
{
  int status = FUNCTION(gsl_block,fread) (stream, m->data, 
					  m->size1 * m->size2,
					  1) ; /* stride */
  return status;
}

int 
FUNCTION(gsl_matrix,fwrite)(FILE * stream, const TYPE(gsl_matrix) * m)
{
  int status = FUNCTION(gsl_block,fwrite) (stream, m->data, 
					   m->size1 * m->size2,
					   1) ; /* stride */
  return status;
}

int
FUNCTION(gsl_matrix,fprintf)(FILE * stream, const TYPE(gsl_matrix) * m,
			     const char * format)
{
  int status = FUNCTION(gsl_block,fprintf) (stream, m->data,
					    m->size1 * m->size2, 
					    1, /* stride */
					    format) ;
  return status;
}

int
FUNCTION(gsl_matrix,fscanf)(FILE * stream, TYPE(gsl_matrix) * m)
{
  int status = FUNCTION(gsl_block,fscanf) (stream, m->data, 
					   m->size1 * m->size2,
					   1) ;  /* stride */
  return status;
}


