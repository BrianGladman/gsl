#ifndef __GSL_BLOCK_COMPLEX_DOUBLE_H__
#define __GSL_BLOCK_COMPLEX_DOUBLE_H__

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

struct gsl_block_complex_struct
{
  size_t size;
  double *data;
};

typedef struct gsl_block_complex_struct gsl_block_complex;

gsl_block_complex *gsl_block_complex_alloc (const size_t n);
gsl_block_complex *gsl_block_complex_calloc (const size_t n);
void gsl_block_complex_free (gsl_block_complex * b);

int gsl_block_complex_fread (FILE * stream, gsl_block_complex * b);
int gsl_block_complex_fwrite (FILE * stream, const gsl_block_complex * b);
int gsl_block_complex_fscanf (FILE * stream, gsl_block_complex * b);
int gsl_block_complex_fprintf (FILE * stream, const gsl_block_complex * b, const char *format);

int gsl_block_complex_raw_fread (FILE * stream, double * b, const size_t n, const size_t stride);
int gsl_block_complex_raw_fwrite (FILE * stream, const double * b, const size_t n, const size_t stride);
int gsl_block_complex_raw_fscanf (FILE * stream, double * b, const size_t n, const size_t stride);
int gsl_block_complex_raw_fprintf (FILE * stream, const double * b, const size_t n, const size_t stride, const char *format);

size_t gsl_block_complex_size (const gsl_block_complex * b);
double * gsl_block_complex_data (const gsl_block_complex * b);

__END_DECLS

#endif /* __GSL_BLOCK_COMPLEX_DOUBLE_H__ */
