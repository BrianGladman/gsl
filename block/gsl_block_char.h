#ifndef GSL_BLOCK_CHAR_H
#define GSL_BLOCK_CHAR_H

#include <stdlib.h>
#include <gsl_errno.h>

struct gsl_block_char_struct
{
  size_t size;
  char *data;
};

typedef struct gsl_block_char_struct gsl_block_char;

gsl_block_char *gsl_block_char_alloc (size_t n);
gsl_block_char *gsl_block_char_calloc (size_t n);
void gsl_block_char_free (gsl_block_char * b);

int gsl_block_char_fread (FILE * stream, gsl_block_char * b);
int gsl_block_char_fwrite (FILE * stream, const gsl_block_char * b);
int gsl_block_char_fscanf (FILE * stream, gsl_block_char * b);
int gsl_block_char_fprintf (FILE * stream, const gsl_block_char * b, const char *format);

int gsl_block_char_fread_with_stride (FILE * stream, char * b, size_t n, size_t stride);
int gsl_block_char_fwrite_with_stride (FILE * stream, const char * b, size_t n, size_t stride);
int gsl_block_char_fscanf_with_stride (FILE * stream, char * b, size_t n, size_t stride);
int gsl_block_char_fprintf_with_stride (FILE * stream, const char * b, size_t n, size_t stride, const char *format);

size_t gsl_block_char_size (const gsl_block_char * b);
char * gsl_block_char_data (const gsl_block_char * b);

#endif /* GSL_BLOCK_CHAR_H */
