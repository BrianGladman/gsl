/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_SF_CHEBYSHEV_H_
#define GSL_SF_CHEBYSHEV_H_


/* data for a Chebyshev series over a given interval */
struct gsl_sf_cheb_series {
  double * c;   /* coefficients */
  int order;    /* order of expansion */
  double a;     /* lower interval point */
  double b;     /* upper interval point */
  double * cp;  /* coefficients of derivative */
  double * ci;  /* coefficients of integral   */
};


/* Calculate a Chebyshev series of specified order over
 * a specified interval, for a given function.
 * Return 0 on failure.
 */
struct gsl_sf_cheb_series * gsl_sf_cheb_new(double (*func)(double),
				            double a, double b,
				            int order
                                            );

/* Calculate a Chebyshev series, but do not allocate
 * a new cheb_series struct. Instead use the one provided.
 * Uses the interval (a,b) and the order with which it
 * was initially created; if you want to change these, then
 * use gsl_sf_cheb_new() instead.
 *
 * exceptions: GSL_EFAULT, GSL_ENOMEM
 */
int gsl_sf_cheb_calc_impl(struct gsl_sf_cheb_series * cs, double (*func)(double));
int gsl_sf_cheb_calc_e(struct gsl_sf_cheb_series * cs, double (*func)(double));


/* Evaluate a Chebyshev series at a given point.
 * No errors can occur for a struct obtained from gsl_sf_cheb_new().
 */
double gsl_sf_cheb_eval(const struct gsl_sf_cheb_series * cs, double x);

/* Evaluate a Chebyshev series at a given point, to (at most) the given order.
 * No errors can occur for a struct obtained from gsl_sf_cheb_new().
 */
double gsl_sf_cheb_eval_n(const struct gsl_sf_cheb_series * cs, int order, double x);

/* Evaluate derivative of a Chebyshev series at a given point.
 */
double gsl_sf_cheb_eval_deriv(struct gsl_sf_cheb_series * cs, double x);

/* Evaluate integal of a Chebyshev series at a given point. The
 * integral is fixed by the condition that it equals zero at
 * the left end-point, ie it is precisely
 *       Integrate[cs(t; a,b), {t, a, x}]
 */
double gsl_sf_cheb_eval_integ(struct gsl_sf_cheb_series * cs, double x);


/* Free a Chebyshev series previously calculated with gsl_sf_cheb_new().
 * No error can occur.
 */
void gsl_sf_cheb_free(struct gsl_sf_cheb_series * cs);


#endif /* GSL_SF_CHEBYSHEV_H_ */
