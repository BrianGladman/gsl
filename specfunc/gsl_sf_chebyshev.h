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
   Return 0 on failure.
 */
struct gsl_sf_ChebSeries * gsl_sf_cheb_new(const double (*func)(double),
				           double a, double b,
				           int order
                                           );

/* Calculate a Chebyshev series, but do not allocate
   a new ChebSeries struct. Instead use the one provided.
   Uses the interval (a,b) and the order with which it
   was initially created; if you want to change these, then
   use gsl_sf_cheb_new() instead.
   Returns a GSL error status, GSL_SUCCESS on success.
 */
int gsl_sf_cheb_calc(struct gsl_sf_ChebSeries * cs, double (*func)(double));


/* Evaluate a Chebyshev series at a given point.
 * No errors can occur for a struct obtained from gsl_sf_cheb_new().
 */
double gsl_sf_cheb_eval(double x, const struct gsl_sf_ChebSeries * cs);


/* Evaluate a Chebyshev series at a given point, to (at most) the given order.
 * No errors can occur for a struct obtained from gsl_sf_cheb_new().
 */
double gsl_sf_cheb_eval_n(double x, int order, const struct gsl_sf_ChebSeries * cs);


/* Free a Chebyshev series previously calculated with gsl_sf_cheb_new().
 * No error can occur.
 */
void gsl_sf_cheb_free(struct gsl_sf_ChebSeries * cs);


#endif /* GSL_SF_CHEBYSHEV_H_ */
