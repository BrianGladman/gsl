#ifndef __GSL_BLOCK_SHORT_H__
#define __GSL_BLOCK_SHORT_H__

#include <stdlib.h>
#include <gsl/gsl_errno.h>

struct gsl_block_short_struct
{
  size_t size;
  short *data;
};

typedef struct gsl_block_short_struct gsl_block_short;

gsl_block_short *gsl_block_short_alloc (size_t n);
gsl_block_short *gsl_block_short_calloc (size_t n);
void gsl_block_short_free (gsl_block_short * b);

int gsl_block_short_fread (FILE * stream, gsl_block_short * b);
int gsl_block_short_fwrite (FILE * stream, const gsl_block_short * b);
int gsl_block_short_fscanf (FILE * stream, gsl_block_short * b);
int gsl_block_short_fprintf (FILE * stream, const gsl_block_short * b, const char *format);

int gsl_block_short_raw_fread (FILE * stream, short * b, size_t n, size_t stride);
int gsl_block_short_raw_fwrite (FILE * stream, const short * b, size_t n, size_t stride);
int gsl_block_short_raw_fscanf (FILE * stream, short * b, size_t n, size_t stride);
int gsl_block_short_raw_fprintf (FILE * stream, const short * b, size_t n, size_t stride, const char *format);

size_t gsl_block_short_size (const gsl_block_short * b);
short * gsl_block_short_data (const gsl_block_short * b);

#endif /* __GSL_BLOCK_SHORT_H__ */
