/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_SF_CHEBYSHEV_H_
#define GSL_SF_CHEBYSHEV_H_

#include <gsl_sf_result.h>


/* data for a Chebyshev series over a given interval */
struct gsl_sf_cheb_series_struct {
  double * c;   /* coefficients                */
  int order;    /* order of expansion          */
  double a;     /* lower interval point        */
  double b;     /* upper interval point        */
  double * cp;  /* coefficients of derivative  */
  double * ci;  /* coefficients of integral    */

  /* the following exists for the benefit of the implementation */
  int order_sp;
};
typedef struct gsl_sf_cheb_series_struct gsl_sf_cheb_series;


/* Calculate a Chebyshev series of specified order over
 * a specified interval, for a given function.
 * Return 0 on failure.
 */
gsl_sf_cheb_series * gsl_sf_cheb_new(double (*func)(double),
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
int gsl_sf_cheb_calc_impl(gsl_sf_cheb_series * cs, double (*func)(double));
int gsl_sf_cheb_calc_e(gsl_sf_cheb_series * cs, double (*func)(double));


/* Evaluate a Chebyshev series at a given point.
 * No errors can occur for a struct obtained from gsl_sf_cheb_new().
 */
int gsl_sf_cheb_eval_impl(const gsl_sf_cheb_series * cs, double x, gsl_sf_result * result);
int gsl_sf_cheb_eval_e(const gsl_sf_cheb_series * cs, double x, gsl_sf_result * result);


/* Evaluate a Chebyshev series at a given point, to (at most) the given order.
 * No errors can occur for a struct obtained from gsl_sf_cheb_new().
 */
int gsl_sf_cheb_eval_n_impl(const gsl_sf_cheb_series * cs, int order, double x, gsl_sf_result * result);
int gsl_sf_cheb_eval_n_e(const gsl_sf_cheb_series * cs, int order, double x, gsl_sf_result * result);


/* Evaluate derivative of a Chebyshev series at a given point.
 */
int gsl_sf_cheb_eval_deriv_impl(gsl_sf_cheb_series * cs, double x, gsl_sf_result * result);
int gsl_sf_cheb_eval_deriv_e(gsl_sf_cheb_series * cs, double x, gsl_sf_result * result);


/* Evaluate integal of a Chebyshev series at a given point. The
 * integral is fixed by the condition that it equals zero at
 * the left end-point, ie it is precisely
 *       Integrate[cs(t; a,b), {t, a, x}]
 */
int gsl_sf_cheb_eval_integ_impl(gsl_sf_cheb_series * cs, double x, gsl_sf_result * result);
int gsl_sf_cheb_eval_integ_e(gsl_sf_cheb_series * cs, double x, gsl_sf_result * result);


/* Free a Chebyshev series previously calculated with gsl_sf_cheb_new().
 * No error can occur.
 */
void gsl_sf_cheb_free(gsl_sf_cheb_series * cs);


#endif /* GSL_SF_CHEBYSHEV_H_ */
