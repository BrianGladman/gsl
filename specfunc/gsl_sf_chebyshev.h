/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_SF_CHEBYSHEV_H_
#define GSL_SF_CHEBYSHEV_H_


/* data for a Chebyshev series over a given interval */
struct gsl_sf_ChebSeries {
  double * c;   /* coefficients */
  int order;    /* order of expansion */
  double a;     /* lower interval point */
  double b;     /* upper interval point */
};


/* Calculate a Chebyshev series of specified order over
   a specified interval, for a given function.
 */
struct gsl_sf_ChebSeries * cheb_new(double (*func)(double),
				    double a, double b,
				    int order);

/* Calculate a Chebyshev series, but do not allocate
   a new ChebSeries. Instead use the one provided.
   The order is assumed to be set to the desired order.
 */
int gsl_sf_cheb_calc(struct gsl_sf_ChebSeries * cs, double (*func)(double));

/* Evaluate a Chebyshev series at a given point. */
double gsl_sf_cheb_eval(double x, const struct gsl_sf_ChebSeries * cs);

/* Evaluate a Chebyshev series at a given point, to (at most) the given order. */
double gsl_sf_cheb_eval_n(double x, int order, const struct gsl_sf_ChebSeries * cs);

/* Free a Chebyshev series previously calculated with gsl_sf_cheb_new(). */
void gsl_sf_cheb_free(struct gsl_sf_ChebSeries * cs);


#endif /* GSL_SF_CHEBYSHEV_H_ */
