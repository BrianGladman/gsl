#ifndef GSL_BLOCK_LONG_H
#define GSL_BLOCK_LONG_H

#include <stdlib.h>
#include <gsl_errno.h>

struct gsl_block_long_struct
{
  size_t size;
  long *data;
};

typedef struct gsl_block_long_struct gsl_block_long;

gsl_block_long *gsl_block_long_alloc (size_t n);
gsl_block_long *gsl_block_long_calloc (size_t n);
void gsl_block_long_free (gsl_block_long * b);

int gsl_block_long_fread (FILE * stream, gsl_block_long * b);
int gsl_block_long_fwrite (FILE * stream, const gsl_block_long * b);
int gsl_block_long_fscanf (FILE * stream, gsl_block_long * b);
int gsl_block_long_fprintf (FILE * stream, const gsl_block_long * b, const char *format);

int gsl_block_long_raw_fread (FILE * stream, long * b, size_t n, size_t stride);
int gsl_block_long_raw_fwrite (FILE * stream, const long * b, size_t n, size_t stride);
int gsl_block_long_raw_fscanf (FILE * stream, long * b, size_t n, size_t stride);
int gsl_block_long_raw_fprintf (FILE * stream, const long * b, size_t n, size_t stride, const char *format);

size_t gsl_block_long_size (const gsl_block_long * b);
long * gsl_block_long_data (const gsl_block_long * b);

#endif /* GSL_BLOCK_LONG_H */
