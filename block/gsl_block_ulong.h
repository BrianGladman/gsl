#ifndef GSL_BLOCK_ULONG_H
#define GSL_BLOCK_ULONG_H

#include <stdlib.h>
#include <gsl_errno.h>

struct gsl_block_ulong_struct
{
  size_t size;
  unsigned long *data;
};

typedef struct gsl_block_ulong_struct gsl_block_ulong;

gsl_block_ulong *gsl_block_ulong_alloc (size_t n);
gsl_block_ulong *gsl_block_ulong_calloc (size_t n);
void gsl_block_ulong_free (gsl_block_ulong * b);

int gsl_block_ulong_fread (FILE * stream, gsl_block_ulong * b);
int gsl_block_ulong_fwrite (FILE * stream, const gsl_block_ulong * b);
int gsl_block_ulong_fscanf (FILE * stream, gsl_block_ulong * b);
int gsl_block_ulong_fprintf (FILE * stream, const gsl_block_ulong * b, const char *format);

int gsl_block_ulong_fread_with_stride (FILE * stream, unsigned long * b, size_t n, size_t stride);
int gsl_block_ulong_fwrite_with_stride (FILE * stream, const unsigned long * b, size_t n, size_t stride);
int gsl_block_ulong_fscanf_with_stride (FILE * stream, unsigned long * b, size_t n, size_t stride);
int gsl_block_ulong_fprintf_with_stride (FILE * stream, const unsigned long * b, size_t n, size_t stride, const char *format);

size_t gsl_block_ulong_size (const gsl_block_ulong * b);
unsigned long * gsl_block_ulong_data (const gsl_block_ulong * b);

#endif /* GSL_BLOCK_ULONG_H */
