#ifndef GSL_BLOCK_UCHAR_H
#define GSL_BLOCK_UCHAR_H

#include <stdlib.h>
#include <gsl_errno.h>

struct gsl_block_uchar_struct
{
  size_t size;
  unsigned char *data;
};

typedef struct gsl_block_uchar_struct gsl_block_uchar;

gsl_block_uchar *gsl_block_uchar_alloc (size_t n);
gsl_block_uchar *gsl_block_uchar_calloc (size_t n);
void gsl_block_uchar_free (gsl_block_uchar * b);

int gsl_block_uchar_fread (FILE * stream, gsl_block_uchar * b);
int gsl_block_uchar_fwrite (FILE * stream, const gsl_block_uchar * b);
int gsl_block_uchar_fscanf (FILE * stream, gsl_block_uchar * b);
int gsl_block_uchar_fprintf (FILE * stream, const gsl_block_uchar * b, const char *format);

int gsl_block_uchar_fread_with_stride (FILE * stream, unsigned char * b, size_t n, size_t stride);
int gsl_block_uchar_fwrite_with_stride (FILE * stream, const unsigned char * b, size_t n, size_t stride);
int gsl_block_uchar_fscanf_with_stride (FILE * stream, unsigned char * b, size_t n, size_t stride);
int gsl_block_uchar_fprintf_with_stride (FILE * stream, const unsigned char * b, size_t n, size_t stride, const char *format);

size_t gsl_block_uchar_size (const gsl_block_uchar * b);
unsigned char * gsl_block_uchar_data (const gsl_block_uchar * b);

#endif /* GSL_BLOCK_UCHAR_H */
