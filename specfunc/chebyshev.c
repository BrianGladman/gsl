#include <stdlib.h>
#include <math.h>
#include <gsl_errno.h>
#include <gsl_math.h>
#include "gsl_sf_chebyshev.h"


struct gsl_sf_ChebSeries * gsl_sf_cheb_new( double (*func)(double),
    	    	    	    	    	    double a, double b,
			      	    	    int order)
{
  int k, j;
  double bma = 0.5 * (b - a);
  double bpa = 0.5 * (b + a);
  double fac = 2./(order +1.);
  double * f = malloc((order+1) * sizeof(double));
  struct gsl_sf_ChebSeries * cs = (struct gsl_sf_ChebSeries *)
    malloc(sizeof(struct gsl_sf_ChebSeries));
  
  if(cs == 0 || f == 0) {
    GSL_ERROR_RETURN("gsl_sf_cheb_new: out of memory",
	    	     GSL_ENOMEM,
	    	     0
	    	     );
  }

  cs->order = order;
  cs->a = a;
  cs->b = b;
  cs->c = malloc((order+1) * sizeof(double));
  if(cs->c == 0) {
    GSL_ERROR_RETURN("gsl_sf_cheb_new: out of memory",
	    	     GSL_ENOMEM,
	    	     0
	    	     );
  }
  
  for(k = 0; k<=order; k++) {
    double y = cos(M_PI * (k+0.5)/(order+1));
    f[k] = func(y*bma + bpa);
  }

  for(j = 0; j<=order; j++) {
    double sum = 0.0;
    for(k = 0; k<=order; k++) sum += f[k]*cos(M_PI * j*(k+0.5)/(order+1));
    cs->c[j] = fac * sum;
  }

  free(f);

  return cs;
}


double gsl_sf_cheb_eval(double x, const struct gsl_sf_ChebSeries * cs)
{
  int j;
  double d  = 0.;
  double dd = 0.;

  double y  = (2.*x - cs->a - cs->b) / (cs->b - cs->a);
  double y2 = 2. * y;

  for(j = cs->order; j>=1; j--) {
    double temp = d;
    d = y2*d - dd + cs->c[j];
    dd = temp;
  }
  return y*d - dd + 0.5 * cs->c[0];
}


void gsl_sf_cheb_free(struct gsl_sf_ChebSeries * cs)
{
  if(cs != 0) {
    if(cs->c != 0) free(cs->c);
    free(cs);
  }
}
