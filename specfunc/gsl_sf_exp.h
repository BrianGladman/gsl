/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_SF_EXP_H_
#define GSL_SF_EXP_H_


/* Provide an exp() function with GSL semantics,
 * i.e. with proper error checking, etc.
 *
 * exceptions: GSL_EOVRFLW, GSL_EUNDRFLW
 */
int     gsl_sf_exp_impl(double x, double * result);
int     gsl_sf_exp_e(double x, double * result);
double  gsl_sf_exp(double x);


/* Exponentiate and apply a given sign:  Exp(x) * Sgn(sgn)
 *
 * exceptions: GSL_EOVRFLW, GSL_EUNDRFLW
 */
int     gsl_sf_exp_sgn_impl(double x, double sgn, double * result);
int     gsl_sf_exp_sgn_e(double x, double sgn, double * result);
double  gsl_sf_exp_sgn(double x, double sgn);


/* Exponentiate and multiply by a given factor:  y * Exp(x)
 *
 * exceptions: GSL_EOVRFLW, GSL_EUNDRFLW
 */
int     gsl_sf_exp_mult_impl(double x, double y, double * result);
int     gsl_sf_exp_mult_e(double x, double y, double * result);
double  gsl_sf_exp_mult(double x, double y);


/* exp(x)-1
 *
 * exceptions: GSL_EOVRFLW
 */
int     gsl_sf_expm1_impl(double x, double * result);
int     gsl_sf_expm1_e(double x, double * result);
double  gsl_sf_expm1(double x);


/* (exp(x)-1)/x = 1 + x/2 + x^2/(2*3) + x^3/(2*3*4) + ...
 *
 * exceptions: GSL_EOVRFLW
 */
int     gsl_sf_exprel_impl(double x, double * result);
int     gsl_sf_exprel_e(double x, double * result);
double  gsl_sf_exprel(double x);


/* 2(exp(x)-1-x)/x^2 = 1 + x/3 + x^2/(3*4) + x^3/(3*4*5) + ...
 *
 * exceptions: GSL_EOVRFLW
 */
int     gsl_sf_exprel_2_impl(double x, double * result);
int     gsl_sf_exprel_2_e(double x, double * result);
double  gsl_sf_exprel_2(double x);


/* Similarly for the N-th generalization of
 * the above. The so-called N-relative exponential
 *
 * exprel_N(x) = N!/x^N (exp(x) - Sum[x^k/k!, {k,0,N-1}])
 *             = 1 + x/(N+1) + x^2/((N+1)(N+2)) + ...
 *             = 1F1(1,1+N,x)
 */
int     gsl_sf_exprel_n_impl(int n, double x, double * result);
int     gsl_sf_exprel_n_e(int n, double x, double * result);
double  gsl_sf_exprel_n(int n, double x);



#ifdef HAVE_INLINE
#include <gsl_math.h>
#include <gsl_errno.h>
extern inline
int gsl_sf_exp_impl(const double x, double * result)
{
  if(x > GSL_LOG_DBL_MAX) {
    *result = 0.0;
    return GSL_EOVRFLW;
  }
  else if(x < GSL_LOG_DBL_MIN) {
    *result = 0.0;
    return GSL_EUNDRFLW;
  }
  else {
    *result = exp(x);
    return GSL_SUCCESS;
  }  
}
extern inline
int gsl_sf_exp_sgn_impl(const double x, const double sgn, double * result)
{
  if(x > GSL_LOG_DBL_MAX) {
    *result = 0.0;
    return GSL_EOVRFLW;
  }
  else if(x < GSL_LOG_DBL_MIN) {
    *result = 0.0;
    return GSL_EUNDRFLW;
  }
  else {
    *result = GSL_SIGN(sgn) * exp(x);
    return GSL_SUCCESS;
  }  
}
#endif  /* HAVE_INLINE */


#endif  /* !GSL_SF_EXP_H_ */
