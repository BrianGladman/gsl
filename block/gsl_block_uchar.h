#ifndef __GSL_BLOCK_UCHAR_H__
#define __GSL_BLOCK_UCHAR_H__

#include <stdlib.h>
#include <gsl/gsl_errno.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

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

int gsl_block_uchar_raw_fread (FILE * stream, unsigned char * b, size_t n, size_t stride);
int gsl_block_uchar_raw_fwrite (FILE * stream, const unsigned char * b, size_t n, size_t stride);
int gsl_block_uchar_raw_fscanf (FILE * stream, unsigned char * b, size_t n, size_t stride);
int gsl_block_uchar_raw_fprintf (FILE * stream, const unsigned char * b, size_t n, size_t stride, const char *format);

size_t gsl_block_uchar_size (const gsl_block_uchar * b);
unsigned char * gsl_block_uchar_data (const gsl_block_uchar * b);

__END_DECLS

#endif /* __GSL_BLOCK_UCHAR_H__ */
