#ifndef GSL_HISTOGRAM_H 
#define GSL_HISTOGRAM_H 

#include <stdlib.h>

typedef struct {
  size_t nbins ;
  double * range ;
  double * bin ;
} gsl_histogram ;

gsl_histogram * gsl_histogram_alloc (size_t n);
gsl_histogram * gsl_histogram_alloc_uniform (size_t n, double xmin, 
					     double xmax);
int gsl_histogram_add (gsl_histogram * h, double x);
int gsl_histogram_accumulate (gsl_histogram * h, double x, double weight);
int gsl_histogram_find (const size_t n, const double * range, 
			const double x, size_t * i);
int gsl_histogram_get (const gsl_histogram * h, size_t i, double * y);

int gsl_histogram_get_binrange (const gsl_histogram * h, size_t i,
				double * x0, double * x1);

double gsl_histogram_max (const gsl_histogram * h);
double gsl_histogram_min (const gsl_histogram * h);
size_t gsl_histogram_nbins (const gsl_histogram * h);

void gsl_histogram_reset (gsl_histogram * h);

#endif /* GSL_HISTOGRAM_H */

