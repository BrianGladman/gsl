#ifndef GSL_BLOCK_UINT_H
#define GSL_BLOCK_UINT_H

#include <stdlib.h>
#include <gsl_errno.h>

struct gsl_block_uint_struct
{
  size_t size;
  unsigned int *data;
};

typedef struct gsl_block_uint_struct gsl_block_uint;

gsl_block_uint *gsl_block_uint_alloc (size_t n);
gsl_block_uint *gsl_block_uint_calloc (size_t n);
void gsl_block_uint_free (gsl_block_uint * b);

int gsl_block_uint_fread (FILE * stream, gsl_block_uint * b);
int gsl_block_uint_fwrite (FILE * stream, const gsl_block_uint * b);
int gsl_block_uint_fscanf (FILE * stream, gsl_block_uint * b);
int gsl_block_uint_fprintf (FILE * stream, const gsl_block_uint * b, const char *format);

int gsl_block_uint_fread_with_stride (FILE * stream, unsigned int * b, size_t n, size_t stride);
int gsl_block_uint_fwrite_with_stride (FILE * stream, const unsigned int * b, size_t n, size_t stride);
int gsl_block_uint_fscanf_with_stride (FILE * stream, unsigned int * b, size_t n, size_t stride);
int gsl_block_uint_fprintf_with_stride (FILE * stream, const unsigned int * b, size_t n, size_t stride, const char *format);

size_t gsl_block_uint_size (const gsl_block_uint * b);
unsigned int * gsl_block_uint_data (const gsl_block_uint * b);

#endif /* GSL_BLOCK_UINT_H */
