#ifndef GSL_MATRIX_CHAR_H
#define GSL_MATRIX_CHAR_H

#include <stdlib.h>
#include <gsl_errno.h>
#include <gsl_config.h>

typedef struct
{
  size_t size1;
  size_t size2;
  char * data;
} gsl_matrix_char ;

gsl_matrix_char * gsl_matrix_char_alloc (size_t n1, size_t n2);
gsl_matrix_char * gsl_matrix_char_calloc (size_t n1, size_t n2);
void gsl_matrix_char_free (gsl_matrix_char * m);

char * gsl_matrix_char_ptr(const gsl_matrix_char * m, size_t i, size_t j);
char   gsl_matrix_char_get(const gsl_matrix_char * m, size_t i, size_t j);
void  gsl_matrix_char_set(gsl_matrix_char * m, size_t i,  size_t j, char x);

int gsl_matrix_char_fread (FILE * stream, gsl_matrix_char * m) ;
int gsl_matrix_char_fwrite (FILE * stream, const gsl_matrix_char * m) ;
int gsl_matrix_char_fscanf (FILE * stream, gsl_matrix_char * m);
int gsl_matrix_char_fprintf (FILE * stream, const gsl_matrix_char * m, const char * format);

extern int gsl_check_range ;


#ifdef HAVE_INLINE

extern inline 
char
gsl_matrix_char_get(const gsl_matrix_char * m, const size_t i, const size_t j)
{
#ifndef GSL_RANGE_CHECK_OFF
  if (i >= m->size1)  /* size_t is unsigned, can't be negative */
    {
      GSL_ERROR_CONTINUE("first index out of range", GSL_EINVAL) ;
    }
  else if (j >= m->size2) /* size_t is unsigned, can't be negative */
    {
      GSL_ERROR_CONTINUE("second index out of range", GSL_EINVAL) ;
    }
#endif
  return m->data[i * m->size2 + j] ;
} 

extern inline 
void
gsl_matrix_char_set(gsl_matrix_char * m, 
		   const size_t i, const size_t j, 
		   const char x)
{
#ifndef GSL_RANGE_CHECK_OFF
  if (i >= m->size1) /* size_t is unsigned, can't be negative */
    {
      GSL_ERROR_RETURN_NOTHING("first index out of range", GSL_EINVAL) ;
    }
  else if (j >= m->size2) /* size_t is unsigned, can't be negative */
    {
      GSL_ERROR_RETURN_NOTHING("second index out of range", GSL_EINVAL) ;
    }
#endif
  m->data[i * m->size2 + j] = x ;
}
#endif /* HAVE_INLINE */

#endif /* !GSL_MATRIX_CHAR_H */
