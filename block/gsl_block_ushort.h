#ifndef __GSL_BLOCK_USHORT_H__
#define __GSL_BLOCK_USHORT_H__

#include <stdlib.h>
#include <gsl/gsl_errno.h>

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

int gsl_block_ushort_raw_fread (FILE * stream, unsigned short * b, size_t n, size_t stride);
int gsl_block_ushort_raw_fwrite (FILE * stream, const unsigned short * b, size_t n, size_t stride);
int gsl_block_ushort_raw_fscanf (FILE * stream, unsigned short * b, size_t n, size_t stride);
int gsl_block_ushort_raw_fprintf (FILE * stream, const unsigned short * b, size_t n, size_t stride, const char *format);

size_t gsl_block_ushort_size (const gsl_block_ushort * b);
unsigned short * gsl_block_ushort_data (const gsl_block_ushort * b);

#endif /* __GSL_BLOCK_USHORT_H__ */
