#ifndef GSL_VECTOR_CHAR_H
#define GSL_VECTOR_CHAR_H

#include <stdlib.h>
#include <gsl_errno.h>
#include <gsl_block_char.h>

struct gsl_vector_char_struct
{
  size_t size;
  size_t stride;
  char *data;
  gsl_block_char *block;
};

typedef struct gsl_vector_char_struct gsl_vector_char;

gsl_vector_char *gsl_vector_char_alloc (size_t n);
gsl_vector_char *gsl_vector_char_calloc (size_t n);

gsl_vector_char *gsl_vector_char_alloc_from_block (gsl_block_char * b,
                                                     size_t offset, 
                                                     size_t n, 
                                                     size_t stride);

gsl_vector_char *gsl_vector_char_alloc_from_vector (gsl_vector_char * v,
                                                      size_t offset, 
                                                      size_t n, 
                                                      size_t stride);

void gsl_vector_char_free (gsl_vector_char * v);

char *gsl_vector_char_ptr (const gsl_vector_char * v, const size_t i);
char gsl_vector_char_get (const gsl_vector_char * v, const size_t i);
void gsl_vector_char_set (gsl_vector_char * v, const size_t i, char x);

int gsl_vector_char_fread (FILE * stream, gsl_vector_char * v);
int gsl_vector_char_fwrite (FILE * stream, const gsl_vector_char * v);
int gsl_vector_char_fscanf (FILE * stream, gsl_vector_char * v);
int gsl_vector_char_fprintf (FILE * stream, const gsl_vector_char * v,
			      const char *format);

int gsl_vector_char_copy (gsl_vector_char * dest, const gsl_vector_char * src);

int gsl_vector_char_swap (gsl_vector_char * v, size_t i, size_t j);

int gsl_vector_char_isnull (gsl_vector_char * v);

extern int gsl_check_range;

#ifdef HAVE_INLINE

extern inline
char
gsl_vector_char_get (const gsl_vector_char * v, const size_t i)
{
#ifndef GSL_RANGE_CHECK_OFF
  if (i >= v->size)
    {
      GSL_ERROR_RETURN ("index out of range", GSL_EINVAL, 0);
    }
#endif
  return v->data[i * v->stride];
}

extern inline
void
gsl_vector_char_set (gsl_vector_char * v, const size_t i, char x)
{
#ifndef GSL_RANGE_CHECK_OFF
  if (i >= v->size)
    {
      GSL_ERROR_RETURN_NOTHING ("index out of range", GSL_EINVAL);
    }
#endif
  v->data[i * v->stride] = x;
}

#endif /* HAVE_INLINE */

#endif /* GSL_VECTOR_CHAR_H */








