#ifndef GSL_BLOCK_COMPLEX_FLOAT_H
#define GSL_BLOCK_COMPLEX_FLOAT_H

#include <stdlib.h>
#include <gsl_errno.h>

struct gsl_block_complex_float_struct
{
  size_t size;
  float *data;
};

typedef struct gsl_block_complex_float_struct gsl_block_complex_float;

gsl_block_complex_float *gsl_block_complex_float_alloc (size_t n);
gsl_block_complex_float *gsl_block_complex_float_calloc (size_t n);
void gsl_block_complex_float_free (gsl_block_complex_float * b);

int gsl_block_complex_float_fread (FILE * stream, gsl_block_complex_float * b);
int gsl_block_complex_float_fwrite (FILE * stream, const gsl_block_complex_float * b);
int gsl_block_complex_float_fscanf (FILE * stream, gsl_block_complex_float * b);
int gsl_block_complex_float_fprintf (FILE * stream, const gsl_block_complex_float * b, const char *format);

int gsl_block_complex_float_fread_with_stride (FILE * stream, float * b, size_t n, size_t stride);
int gsl_block_complex_float_fwrite_with_stride (FILE * stream, const float * b, size_t n, size_t stride);
int gsl_block_complex_float_fscanf_with_stride (FILE * stream, float * b, size_t n, size_t stride);
int gsl_block_complex_float_fprintf_with_stride (FILE * stream, const float * b, size_t n, size_t stride, const char *format);

size_t gsl_block_complex_float_size (const gsl_block_complex_float * b);
float * gsl_block_complex_float_data (const gsl_block_complex_float * b);

#endif /* GSL_BLOCK_COMPLEX_FLOAT_H */
