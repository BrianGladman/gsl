#ifndef __GSL_BLOCK_DOUBLE_H__
#define __GSL_BLOCK_DOUBLE_H__

#include <stdlib.h>
#include <gsl/gsl_errno.h>

struct gsl_block_struct
{
  size_t size;
  double *data;
};

typedef struct gsl_block_struct gsl_block;

gsl_block *gsl_block_alloc (size_t n);
gsl_block *gsl_block_calloc (size_t n);
void gsl_block_free (gsl_block * b);

int gsl_block_fread (FILE * stream, gsl_block * b);
int gsl_block_fwrite (FILE * stream, const gsl_block * b);
int gsl_block_fscanf (FILE * stream, gsl_block * b);
int gsl_block_fprintf (FILE * stream, const gsl_block * b, const char *format);

int gsl_block_raw_fread (FILE * stream, double * b, size_t n, size_t stride);
int gsl_block_raw_fwrite (FILE * stream, const double * b, size_t n, size_t stride);
int gsl_block_raw_fscanf (FILE * stream, double * b, size_t n, size_t stride);
int gsl_block_raw_fprintf (FILE * stream, const double * b, size_t n, size_t stride, const char *format);

size_t gsl_block_size (const gsl_block * b);
double * gsl_block_data (const gsl_block * b);

#endif /* __GSL_BLOCK_DOUBLE_H__ */
