#ifndef GSL_BLOCK_COMPLEX_DOUBLE_H
#define GSL_BLOCK_COMPLEX_DOUBLE_H

#include <stdlib.h>
#include <gsl_errno.h>

struct gsl_block_complex_struct
{
  size_t size;
  double *data;
};

typedef struct gsl_block_complex_struct gsl_block_complex;

gsl_block_complex *gsl_block_complex_alloc (size_t n);
gsl_block_complex *gsl_block_complex_calloc (size_t n);
void gsl_block_complex_free (gsl_block_complex * b);

int gsl_block_complex_fread (FILE * stream, gsl_block_complex * b);
int gsl_block_complex_fwrite (FILE * stream, const gsl_block_complex * b);
int gsl_block_complex_fscanf (FILE * stream, gsl_block_complex * b);
int gsl_block_complex_fprintf (FILE * stream, const gsl_block_complex * b, const char *format);

int gsl_block_complex_raw_fread (FILE * stream, double * b, size_t n, size_t stride);
int gsl_block_complex_raw_fwrite (FILE * stream, const double * b, size_t n, size_t stride);
int gsl_block_complex_raw_fscanf (FILE * stream, double * b, size_t n, size_t stride);
int gsl_block_complex_raw_fprintf (FILE * stream, const double * b, size_t n, size_t stride, const char *format);

size_t gsl_block_complex_size (const gsl_block_complex * b);
double * gsl_block_complex_data (const gsl_block_complex * b);

#endif /* GSL_BLOCK_COMPLEX_DOUBLE_H */
