#ifndef __GSL_HISTOGRAM_H__
#define __GSL_HISTOGRAM_H__

#include <stdlib.h>
#include <stdio.h>

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

typedef struct {
  size_t n ;
  double * range ;
  double * bin ;
} gsl_histogram ;

typedef struct {
  size_t n ;
  double * range ;
  double * sum ;
} gsl_histogram_pdf ;

gsl_histogram * gsl_histogram_calloc (size_t n);
gsl_histogram * gsl_histogram_calloc_uniform (const size_t n, const double xmin, 
					     const double xmax);
void gsl_histogram_free (gsl_histogram * h);
int gsl_histogram_increment (gsl_histogram * h, double x);
int gsl_histogram_accumulate (gsl_histogram * h, double x, double weight);
int gsl_histogram_find (const gsl_histogram * h, 
			const double x, size_t * i);
int gsl_histogram_find_impl (const size_t n, const double range[],
			     const double x, size_t * i);

double gsl_histogram_get (const gsl_histogram * h, size_t i);
int gsl_histogram_get_range (const gsl_histogram * h, size_t i, 
			     double * lower, double * upper);
				     
double gsl_histogram_max (const gsl_histogram * h);
double gsl_histogram_min (const gsl_histogram * h);
size_t gsl_histogram_bins (const gsl_histogram * h);

void gsl_histogram_reset (gsl_histogram * h);

int gsl_histogram_fwrite (FILE * stream, const gsl_histogram * h) ;
int gsl_histogram_fread (FILE * stream, gsl_histogram * h);
int gsl_histogram_fprintf (FILE * stream, const gsl_histogram * h, 
			   const char * range_format, const char * bin_format);
int gsl_histogram_fscanf (FILE * stream, gsl_histogram * h);

gsl_histogram_pdf * gsl_histogram_pdf_alloc (const gsl_histogram * h);
void gsl_histogram_pdf_free (gsl_histogram_pdf * p);
double gsl_histogram_pdf_sample (const gsl_histogram_pdf * p, double r);

__END_DECLS

#endif /* __GSL_HISTOGRAM_H__ */
