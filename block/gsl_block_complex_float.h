#ifndef __GSL_BLOCK_COMPLEX_FLOAT_H__
#define __GSL_BLOCK_COMPLEX_FLOAT_H__

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

struct gsl_block_complex_float_struct
{
  size_t size;
  float *data;
};

typedef struct gsl_block_complex_float_struct gsl_block_complex_float;

gsl_block_complex_float *gsl_block_complex_float_alloc (const size_t n);
gsl_block_complex_float *gsl_block_complex_float_calloc (const size_t n);
void gsl_block_complex_float_free (gsl_block_complex_float * b);

int gsl_block_complex_float_fread (FILE * stream, gsl_block_complex_float * b);
int gsl_block_complex_float_fwrite (FILE * stream, const gsl_block_complex_float * b);
int gsl_block_complex_float_fscanf (FILE * stream, gsl_block_complex_float * b);
int gsl_block_complex_float_fprintf (FILE * stream, const gsl_block_complex_float * b, const char *format);

int gsl_block_complex_float_raw_fread (FILE * stream, float * b, const size_t n, const size_t stride);
int gsl_block_complex_float_raw_fwrite (FILE * stream, const float * b, const size_t n, const size_t stride);
int gsl_block_complex_float_raw_fscanf (FILE * stream, float * b, const size_t n, const size_t stride);
int gsl_block_complex_float_raw_fprintf (FILE * stream, const float * b, const size_t n, const size_t stride, const char *format);

size_t gsl_block_complex_float_size (const gsl_block_complex_float * b);
float * gsl_block_complex_float_data (const gsl_block_complex_float * b);

__END_DECLS

#endif /* __GSL_BLOCK_COMPLEX_FLOAT_H__ */
