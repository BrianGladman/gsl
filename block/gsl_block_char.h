#ifndef __GSL_BLOCK_CHAR_H__
#define __GSL_BLOCK_CHAR_H__

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

int gsl_block_char_raw_fread (FILE * stream, char * b, size_t n, size_t stride);
int gsl_block_char_raw_fwrite (FILE * stream, const char * b, size_t n, size_t stride);
int gsl_block_char_raw_fscanf (FILE * stream, char * b, size_t n, size_t stride);
int gsl_block_char_raw_fprintf (FILE * stream, const char * b, size_t n, size_t stride, const char *format);

size_t gsl_block_char_size (const gsl_block_char * b);
char * gsl_block_char_data (const gsl_block_char * b);

__END_DECLS

#endif /* __GSL_BLOCK_CHAR_H__ */
