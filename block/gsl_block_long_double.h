#ifndef __GSL_BLOCK_LONG_DOUBLE_H__
#define __GSL_BLOCK_LONG_DOUBLE_H__

#include <stdlib.h>
#include <gsl/gsl_errno.h>

struct gsl_block_long_double_struct
{
  size_t size;
  long double *data;
};

typedef struct gsl_block_long_double_struct gsl_block_long_double;

gsl_block_long_double *gsl_block_long_double_alloc (size_t n);
gsl_block_long_double *gsl_block_long_double_calloc (size_t n);
void gsl_block_long_double_free (gsl_block_long_double * b);

int gsl_block_long_double_fread (FILE * stream, gsl_block_long_double * b);
int gsl_block_long_double_fwrite (FILE * stream, const gsl_block_long_double * b);
int gsl_block_long_double_fscanf (FILE * stream, gsl_block_long_double * b);
int gsl_block_long_double_fprintf (FILE * stream, const gsl_block_long_double * b, const char *format);

int gsl_block_long_double_raw_fread (FILE * stream, long double * b, size_t n, size_t stride);
int gsl_block_long_double_raw_fwrite (FILE * stream, const long double * b, size_t n, size_t stride);
int gsl_block_long_double_raw_fscanf (FILE * stream, long double * b, size_t n, size_t stride);
int gsl_block_long_double_raw_fprintf (FILE * stream, const long double * b, size_t n, size_t stride, const char *format);

size_t gsl_block_long_double_size (const gsl_block_long_double * b);
long double * gsl_block_long_double_data (const gsl_block_long_double * b);

#endif /* __GSL_BLOCK_LONG_DOUBLE_H__ */
