/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef __GSL_SF_CHEBYSHEV_H__
#define __GSL_SF_CHEBYSHEV_H__

#include <gsl/gsl_mode.h>
#include <gsl/gsl_sf_result.h>

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


/* data for a Chebyshev series over a given interval */
struct gsl_sf_cheb_series_struct {
  double * c;   /* coefficients                */
  int order;    /* order of expansion          */
  double a;     /* lower interval point        */
  double b;     /* upper interval point        */
  double * cp;  /* coefficients of derivative  */
  double * ci;  /* coefficients of integral    */

  /* The following exists (mostly) for the benefit
   * of the implementation. It is an effective single
   * precision order, for use in single precision
   * evaluation. Users can use it if they like, but
   * only they know how to calculate it, since it is
   * specific to the approximated function. By default,
   * order_sp = order.
   * It is used explicitly only by the gsl_sf_cheb_eval_mode
   * functions, which are not meant for casual use.
   */
  int order_sp;
};
typedef struct gsl_sf_cheb_series_struct gsl_sf_cheb_series;


/* Calculate a Chebyshev series of specified order over
 * a specified interval, for a given function.
 * Return 0 on failure.
 */
gsl_sf_cheb_series * gsl_sf_cheb_new(double (*func)(double),
				     const double a, const double b,
				     const int order
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
int gsl_sf_cheb_eval_impl(const gsl_sf_cheb_series * cs, const double x, gsl_sf_result * result);
int gsl_sf_cheb_eval_e(const gsl_sf_cheb_series * cs, double x, gsl_sf_result * result);


/* Evaluate a Chebyshev series at a given point, to (at most) the given order.
 * No errors can occur for a struct obtained from gsl_sf_cheb_new().
 */
int gsl_sf_cheb_eval_n_impl(const gsl_sf_cheb_series * cs, const int order, const double x, gsl_sf_result * result);
int gsl_sf_cheb_eval_n_e(const gsl_sf_cheb_series * cs, int order, double x, gsl_sf_result * result);


/* Evaluate a Chebyshev series at a given point, using the default
 * order for double precision mode(s) and the single precision
 * order for other modes.
 * No errors can occur for a struct obtained from gsl_sf_cheb_new().
 */
int gsl_sf_cheb_eval_mode_impl(const gsl_sf_cheb_series * cs, const double x, gsl_mode_t mode, gsl_sf_result * result);
int gsl_sf_cheb_eval_mode_e(const gsl_sf_cheb_series * cs, double x, gsl_mode_t mode, gsl_sf_result * result);


/* Evaluate derivative of a Chebyshev series at a given point.
 */
int gsl_sf_cheb_eval_deriv_impl(gsl_sf_cheb_series * cs, const double x, gsl_sf_result * result);
int gsl_sf_cheb_eval_deriv_e(gsl_sf_cheb_series * cs, double x, gsl_sf_result * result);


/* Evaluate integal of a Chebyshev series at a given point. The
 * integral is fixed by the condition that it equals zero at
 * the left end-point, ie it is precisely
 *       Integrate[cs(t; a,b), {t, a, x}]
 */
int gsl_sf_cheb_eval_integ_impl(gsl_sf_cheb_series * cs, const double x, gsl_sf_result * result);
int gsl_sf_cheb_eval_integ_e(gsl_sf_cheb_series * cs, double x, gsl_sf_result * result);


/* Free a Chebyshev series previously calculated with gsl_sf_cheb_new().
 * No error can occur.
 */
void gsl_sf_cheb_free(gsl_sf_cheb_series * cs);


__END_DECLS

#endif /* __GSL_SF_CHEBYSHEV_H__ */
