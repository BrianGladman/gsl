#ifndef __GSL_BLOCK_COMPLEX_LONG_DOUBLE_H__
#define __GSL_BLOCK_COMPLEX_LONG_DOUBLE_H__

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

struct gsl_block_complex_long_double_struct
{
  size_t size;
  long double *data;
};

typedef struct gsl_block_complex_long_double_struct gsl_block_complex_long_double;

gsl_block_complex_long_double *gsl_block_complex_long_double_alloc (const size_t n);
gsl_block_complex_long_double *gsl_block_complex_long_double_calloc (const size_t n);
void gsl_block_complex_long_double_free (gsl_block_complex_long_double * b);

int gsl_block_complex_long_double_fread (FILE * stream, gsl_block_complex_long_double * b);
int gsl_block_complex_long_double_fwrite (FILE * stream, const gsl_block_complex_long_double * b);
int gsl_block_complex_long_double_fscanf (FILE * stream, gsl_block_complex_long_double * b);
int gsl_block_complex_long_double_fprintf (FILE * stream, const gsl_block_complex_long_double * b, const char *format);

int gsl_block_complex_long_double_raw_fread (FILE * stream, long double * b, const size_t n, const size_t stride);
int gsl_block_complex_long_double_raw_fwrite (FILE * stream, const long double * b, const size_t n, const size_t stride);
int gsl_block_complex_long_double_raw_fscanf (FILE * stream, long double * b, const size_t n, const size_t stride);
int gsl_block_complex_long_double_raw_fprintf (FILE * stream, const long double * b, const size_t n, const size_t stride, const char *format);

size_t gsl_block_complex_long_double_size (const gsl_block_complex_long_double * b);
long double * gsl_block_complex_long_double_data (const gsl_block_complex_long_double * b);

__END_DECLS

#endif /* __GSL_BLOCK_COMPLEX_LONG_DOUBLE_H__ */
