/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_SF_EXP_H_
#define GSL_SF_EXP_H_

#include <gsl_sf_result.h>
#include <gsl_precision.h>


/* Provide an exp() function with GSL semantics,
 * i.e. with proper error checking, etc.
 *
 * exceptions: GSL_EOVRFLW, GSL_EUNDRFLW
 */
int gsl_sf_exp_impl(double x, gsl_sf_result * result);
int gsl_sf_exp_e(double x, gsl_sf_result * result);


/* Exponentiate and apply a given sign:  Exp(x) * Sgn(sgn)
 *
 * exceptions: GSL_EOVRFLW, GSL_EUNDRFLW
 */
int gsl_sf_exp_sgn_impl(double x, double sgn, gsl_sf_result * result);
int gsl_sf_exp_sgn_e(double x, double sgn, gsl_sf_result * result);


/* Exponentiate and multiply by a given factor:  y * Exp(x)
 *
 * exceptions: GSL_EOVRFLW, GSL_EUNDRFLW
 */
int gsl_sf_exp_mult_impl(double x, double y, gsl_sf_result * result);
int gsl_sf_exp_mult_e(double x, double y, gsl_sf_result * result);


/* exp(x)-1
 *
 * exceptions: GSL_EOVRFLW
 */
int gsl_sf_expm1_impl(double x, gsl_sf_result * result);
int gsl_sf_expm1_e(double x, gsl_sf_result * result);


/* (exp(x)-1)/x = 1 + x/2 + x^2/(2*3) + x^3/(2*3*4) + ...
 *
 * exceptions: GSL_EOVRFLW
 */
int gsl_sf_exprel_impl(double x, gsl_sf_result * result);
int gsl_sf_exprel_e(double x, gsl_sf_result * result);


/* 2(exp(x)-1-x)/x^2 = 1 + x/3 + x^2/(3*4) + x^3/(3*4*5) + ...
 *
 * exceptions: GSL_EOVRFLW
 */
int gsl_sf_exprel_2_impl(double x, gsl_sf_result * result);
int gsl_sf_exprel_2_e(double x, gsl_sf_result * result);


/* Similarly for the N-th generalization of
 * the above. The so-called N-relative exponential
 *
 * exprel_N(x) = N!/x^N (exp(x) - Sum[x^k/k!, {k,0,N-1}])
 *             = 1 + x/(N+1) + x^2/((N+1)(N+2)) + ...
 *             = 1F1(1,1+N,x)
 */
int gsl_sf_exprel_n_impl(int n, double x, gsl_sf_result * result);
int gsl_sf_exprel_n_e(int n, double x, gsl_sf_result * result);


/* Exponentiate a quantity with an associated error.
 */
int gsl_sf_exp_err_impl(double x, double dx, gsl_sf_result * result);



#ifdef HAVE_INLINE
#include <gsl_math.h>
#include <gsl_errno.h>
extern inline
int gsl_sf_exp_impl(const double x, gsl_sf_result * result)
{
  if(x > GSL_LOG_DBL_MAX) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_EOVRFLW;
  }
  else if(x < GSL_LOG_DBL_MIN) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_EUNDRFLW;
  }
  else {
    result->val = exp(x);
    result->err = GSL_DBL_EPSILON * result->val;
    return GSL_SUCCESS;
  }  
}
extern inline
int gsl_sf_exp_sgn_impl(const double x, const double sgn, gsl_sf_result * result)
{
  if(x > GSL_LOG_DBL_MAX) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_EOVRFLW;
  }
  else if(x < GSL_LOG_DBL_MIN) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_EUNDRFLW;
  }
  else {
    double s = GSL_SIGN(sgn);
    double e = exp(x);
    result->val = s * e;
    result->err = e * GSL_DBL_EPSILON;
    return GSL_SUCCESS;
  }  
}
#endif  /* HAVE_INLINE */


#endif  /* !GSL_SF_EXP_H_ */
