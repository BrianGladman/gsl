#ifndef GSL_BLOCK_INT_H
#define GSL_BLOCK_INT_H

#include <stdlib.h>
#include <gsl/gsl_errno.h>

struct gsl_block_int_struct
{
  size_t size;
  int *data;
};

typedef struct gsl_block_int_struct gsl_block_int;

gsl_block_int *gsl_block_int_alloc (size_t n);
gsl_block_int *gsl_block_int_calloc (size_t n);
void gsl_block_int_free (gsl_block_int * b);

int gsl_block_int_fread (FILE * stream, gsl_block_int * b);
int gsl_block_int_fwrite (FILE * stream, const gsl_block_int * b);
int gsl_block_int_fscanf (FILE * stream, gsl_block_int * b);
int gsl_block_int_fprintf (FILE * stream, const gsl_block_int * b, const char *format);

int gsl_block_int_raw_fread (FILE * stream, int * b, size_t n, size_t stride);
int gsl_block_int_raw_fwrite (FILE * stream, const int * b, size_t n, size_t stride);
int gsl_block_int_raw_fscanf (FILE * stream, int * b, size_t n, size_t stride);
int gsl_block_int_raw_fprintf (FILE * stream, const int * b, size_t n, size_t stride, const char *format);

size_t gsl_block_int_size (const gsl_block_int * b);
int * gsl_block_int_data (const gsl_block_int * b);

#endif /* GSL_BLOCK_INT_H */
