#ifndef GSL_BLOCK_USHORT_H
#define GSL_BLOCK_USHORT_H

#include <stdlib.h>
#include <gsl_errno.h>

struct gsl_block_ushort_struct
{
  size_t size;
  unsigned short *data;
};

typedef struct gsl_block_ushort_struct gsl_block_ushort;

gsl_block_ushort *gsl_block_ushort_alloc (size_t n);
gsl_block_ushort *gsl_block_ushort_calloc (size_t n);
void gsl_block_ushort_free (gsl_block_ushort * b);

int gsl_block_ushort_fread (FILE * stream, gsl_block_ushort * b);
int gsl_block_ushort_fwrite (FILE * stream, const gsl_block_ushort * b);
int gsl_block_ushort_fscanf (FILE * stream, gsl_block_ushort * b);
int gsl_block_ushort_fprintf (FILE * stream, const gsl_block_ushort * b, const char *format);

int gsl_block_ushort_fread_with_stride (FILE * stream, unsigned short * b, size_t n, size_t stride);
int gsl_block_ushort_fwrite_with_stride (FILE * stream, const unsigned short * b, size_t n, size_t stride);
int gsl_block_ushort_fscanf_with_stride (FILE * stream, unsigned short * b, size_t n, size_t stride);
int gsl_block_ushort_fprintf_with_stride (FILE * stream, const unsigned short * b, size_t n, size_t stride, const char *format);

size_t gsl_block_ushort_size (const gsl_block_ushort * b);
unsigned short * gsl_block_ushort_data (const gsl_block_ushort * b);

#endif /* GSL_BLOCK_USHORT_H */
