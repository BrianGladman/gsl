#ifndef GSL_MATRIX_UCHAR_H
#define GSL_MATRIX_UCHAR_H

#include <stdlib.h>
#include <gsl_errno.h>
#include <gsl_config.h>
#include <gsl_vector_uchar.h>

typedef struct gsl_matrix_uchar_struct gsl_matrix_uchar;

struct gsl_matrix_uchar_struct
{
  size_t size1;
  size_t size2;
  unsigned char * data;
} ;

gsl_matrix_uchar * gsl_matrix_uchar_alloc (size_t n1, size_t n2);
gsl_matrix_uchar * gsl_matrix_uchar_calloc (size_t n1, size_t n2);
void gsl_matrix_uchar_free (gsl_matrix_uchar * m);

unsigned char * gsl_matrix_uchar_ptr(const gsl_matrix_uchar * m, size_t i, size_t j);
unsigned char   gsl_matrix_uchar_get(const gsl_matrix_uchar * m, size_t i, size_t j);
void  gsl_matrix_uchar_set(gsl_matrix_uchar * m, size_t i,  size_t j, unsigned char x);

int gsl_matrix_uchar_fread (FILE * stream, gsl_matrix_uchar * m) ;
int gsl_matrix_uchar_fwrite (FILE * stream, const gsl_matrix_uchar * m) ;
int gsl_matrix_uchar_fscanf (FILE * stream, gsl_matrix_uchar * m);
int gsl_matrix_uchar_fprintf (FILE * stream, const gsl_matrix_uchar * m, const char * format);

int gsl_matrix_uchar_copy_row(const gsl_matrix_uchar * m, size_t i, gsl_vector_uchar * v);
int gsl_matrix_uchar_copy_col(const gsl_matrix_uchar * m, size_t j, gsl_vector_uchar * v);
int gsl_matrix_uchar_set_row(gsl_matrix_uchar * m, size_t i, const gsl_vector_uchar * v);
int gsl_matrix_uchar_set_col(gsl_matrix_uchar * m, size_t j, const gsl_vector_uchar * v);

extern int gsl_check_range ;

#ifdef HAVE_INLINE

extern inline 
unsigned char
gsl_matrix_uchar_get(const gsl_matrix_uchar * m, const size_t i, const size_t j)
{
#ifndef GSL_RANGE_CHECK_OFF
  if (i >= m->size1)  /* size_t is unsigned, can't be negative */
    {
      GSL_ERROR_RETURN("first index out of range", GSL_EINVAL, 0) ;
    }
  else if (j >= m->size2) /* size_t is unsigned, can't be negative */
    {
      GSL_ERROR_RETURN("second index out of range", GSL_EINVAL, 0) ;
    }
#endif
  return m->data[i * m->size2 + j] ;
} 

extern inline 
void
gsl_matrix_uchar_set(gsl_matrix_uchar * m, 
		   const size_t i, const size_t j, 
		   const unsigned char x)
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

#endif /* !GSL_MATRIX_UCHAR_H */
