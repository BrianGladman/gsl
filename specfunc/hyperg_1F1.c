/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "hyperg.h"
#include "gsl_sf_exp.h"
#include "gsl_sf_hyperg.h"

#define locMAX(a,b)     ((a) > (b) ? (a) : (b))
#define locEPS          (1000.0*GSL_MACH_EPS)



/* Do the (stable) upward recursion on the parameter 'a'.
 * Work in terms of the function
 *     Y(n) := Gamma(a+n)/Gamma(1+a+n-b) 1F1(a+n;b;x)
 *
 *     (a+n+1-b)Y(n+1) + (b-2a-2n-x)Y(n) + (a+n-1)Y(n-1) = 0
 *
 * Note that this will bomb if a-b+1 is a negative integer and
 * the recursion passes through that integer.
 */
static
int
hyperg_1F1_Y_recurse_a(double a, double b, double x,
                       double n, double Ynm1, double Yn,
		       double N,
		       double * YN)
{
  int k;
  double Ykm1 = Ynm1;
  double Yk   = Yn;
  double Ykp1;

  for(k=n; k<N; k++) {
    double ckp1 = (a + k + 1 - b);
    double ck   = b - 2*a - 2*k - x;
    double ckm1 = a + k - 1;
    
    if(fabs(ckp1) < GSL_MACH_EPS) {
      *YN = 0.0;
      return GSL_EDOM;
    }
    Ykp1 = (-ck*Yk - ckm1*Ykm1)/ckp1;
    
    Ykm1 = Yk;
    Yk   = Ykp1;
  }
  
  *YN = Ykp1;
  return GSL_SUCCESS;
}


/* Manage the upward recursion on the parameter 'a',
 * Evaluating 1F1(a+n_stop,b,x) by recursing up
 * from 1F1(a+n_start,b,x) and 1F1(a+n_start+1,b,x).
 *
 * Assumes n_stop > n_start.
 *
 * Uses the series representation to evaluate at the 
 * reduced values of 'a'. This can be very inefficient
 * if x is large. So it is better not to use this for
 * large x.
 */
static
int
hyperg_1F1_recurse_a(double a, double b, double x,
                     int n_start,
                     int n_stop,
                     double * result)
{
  double prec;
  double Y0, Y1;
  double F0, F1;

  int stat_0 = gsl_sf_hyperg_1F1_series_impl(a+n_start,     b, x, &F0, &prec);
  int stat_1 = gsl_sf_hyperg_1F1_series_impl(a+n_start+1.0, b, x, &F1, &prec);
  
  double lg_a0;       /* log(Gamma(a+n_start))     */
  double lg_ab0;      /* log(Gamma(1+a+n_start-b)) */
  int stat_lg_a0  = gsl_sf_lngamma_impl(a+n_start,     &lg_a0);
  int stat_lg_ab0 = gsl_sf_lngamma_impl(1+a+n_start-b, &lg_ab0);

  


  
  
}


static
int
hyperg_1F1_asymp(const double a, const double b, const double x,
                 double * result, double * prec
                 )
{
  double pre, ln_pre, F;

  if(x == 0.0) {
    *prec   = 1.0;
    *result = 0.0; /* FIXME: ?? */
    return GSL_ELOSS;
  }
  else if(x < 0.0) {
    ln_pre = gsl_sf_lngamma(b) - a*log(-x) - gsl_sf_lngamma(b-a);
  }
  else {
    ln_pre = gsl_sf_lngamma(b) + x + (a-b)*log(x) - gsl_sf_lngamma(a);
  }
  
  
}


int
gsl_sf_hyperg_1F1_impl(const double a, const double b, const double x,
                       double * result
                       )
{
  int a_neg_integer;    /*  a   negative integer  */
  int b_neg_integer;    /*  b   negative integer  */
  int bma_neg_integer;  /*  b-a negative integer  */
  int amb_neg_integer;  /*  a-b negative integer  */

  const double bma = b - a;
  const double amb = a - b;

  a_neg_integer = ( a < 0.0  &&  fabs(a - rint(a)) < locEPS );
  b_neg_integer = ( b < 0.0  &&  fabs(b - rint(b)) < locEPS );
  bma_neg_integer = ( bma < 0.0  &&  fabs(bma - rint(bma)) < locEPS );
  amb_neg_integer = ( amb < 0.0  &&  fabs(amb - rint(amb)) < locEPS );
  
  /* case: a==b,  exp(x) */
  if(fabs(b-a) < locEPS) {
    return gsl_sf_exp_impl(x, result);
  }

  /* case: denominator zeroes before numerator */
  if(b_neg_integer && !(a_neg_integer && a > b + 0.1)) {
    *result = 0.0;
    return GSL_EDOM;
  }

  /* If 'a' is a negative integer, then the
   * series truncates to a polynomial.
   * We do not have to worry about negative integer 'b',
   * since that error condition is trapped above.
   */
  if(a_neg_integer) {
    double prec;
    return gsl_sf_hyperg_1F1_series_impl(a, b, x, result, &prec);
  }

  /* If b-a is a negative integer, use the Kummer transformation
   *    1F1(a,b,x) = Exp(x) 1F1(b-a,b,x)
   * to reduce it to a polynomial times an exponential.
   * Note that there can be no error condition here, since
   * there can only be an error if 'b' is a negative integer, but
   * in that case we would not have gotten this far unless 'a' was
   * a negative integer as well, in which case the above code block
   * handled the situation.
   */
  if(bma_neg_integer) {
    double prec;
    double Ex, Kummer_1F1;
    double int_bma = floor(bma + 0.1);
    int stat_E = gsl_sf_exp_impl(x, &Ex);
    int stat_K = gsl_sf_hyperg_1F1_series_impl(int_bma, b, -x, &Kummer_1F1, &prec);
    int stat;
    if(stat_E != GSL_SUCCESS) {
      *result = 0.0;
      stat    = stat_E;
    }
    else if(stat_K == GSL_ELOSS || prec > locEPS) {
      *result = Ex * Kummer_1F1;
      stat    = GSL_ELOSS;
    }
    else if(stat_K != GSL_SUCCESS) {
      *result = 0.0;
      stat    = stat_K;
    }
    else {
      *result = Ex * Kummer_1F1;
      stat = GSL_SUCCESS;
    }
    return stat;
  }


  /* Now we have dealt with any special negative integer cases,
   * including the error cases, so we are left with a well-defined
   * series evaluation, though the arguments may be large.
   */

  if(fabs(x) < 20.0) {
    if(b > fabs(a) || fabs(a) < 20.0) {
      double prec;
      return gsl_sf_hyperg_1F1_series_impl(a, b, x, result, &prec);
    }
    else if(b > fabs(b-a) || fabs(b-a) < 20.0) {
      double prec;
      double Ex = exp(x);
      double Kummer_1F1;
      int stat_Kummer = gsl_sf_hyperg_1F1_series_impl(b-a, b, x, &Kummer_1F1, &prec);
      *result = Ex * Kummer_1F1;
      return stat_Kummer;
    }
    else {
      /* |a| >> |b|  && (|a| > 30 || |b| > 30) */
      
      if(a < 0.0) {
        /* apply Kummer transformation to */
      }
      else {
      }
    }
  }
  else {
  }

}


int
gsl_sf_hyperg_1F1_e(double a, double b, double x, double * result)
{
  int status = gsl_sf_hyperg_1F1_impl(a, b, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_hyperg_1F1_e", status);
  }
  return status;
}

double
gsl_sf_hyperg_1F1(double a, double b, double x)
{
  double y;
  int status = gsl_sf_hyperg_1F1_impl(a, b, x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_hyperg_1F1", status);
  }
  return y;
}
