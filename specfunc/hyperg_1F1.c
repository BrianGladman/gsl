/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_exp.h"
#include "gsl_sf_hyperg.h"

#define locMAX(a,b)     ((a) > (b) ? (a) : (b))
#define locEPS          (1000.0*GSL_MACH_EPS)



static int hyperg_1F1_series(const double a, const double b, const double x,
                             double * result, double * prec
			     )
{
  double an  = a;
  double bn  = b;
  double n   = 1.0;
  double sum = 1.0;
  double del = 1.0;
  double abs_del = 1.0;
  double max_abs_del = 1.0;
  
  while(abs_del > GSL_MACH_EPS) {
    double u, abs_u;

    if(bn == 0.0) {
      return GSL_EDOM;
    }
    if(an == 0.0 || n > 200.0) {
      max_abs_del *= GSL_MACH_EPS;
      *prec   = fabs(GSL_MACH_EPS * n + max_abs_del);
      *result = sum;
      return GSL_SUCCESS;
    }
    
    u = x * (an/(bn*n));
    abs_u = fabs(u);
    if(abs_u > 1.0 && max_abs_del > DBL_MAX/abs_u) {
      *prec   = 1.0;
      *result = sum;
      return GSL_ELOSS;
    }
    del *= u;
    sum += del;

    max_abs_del = locMAX(fabs(del), max_abs_del);

    an += 1.0;
    bn += 1.0;
    n  += 1.0;
  }
  
  max_abs_del *= GSL_MACH_EPS;
  *prec   = fabs(GSL_MACH_EPS * n + max_abs_del);
  *result = sum;
  if(*prec > locEPS) {
    return GSL_ELOSS;
  }
  else {
    return GSL_SUCCESS;
  }
}

/* [Carlson, p.109] says the error in truncating this asymptotic series
 * is less than the absolute value of the first neglected term
 */
static int hyperg_2F0_series(const double a, const double b, const double x,
                             double * result, double * prec
			     )
{
  double an = a;
  double bn = b;
}

static int hyperg_1F1_asymp(const double a, const double b, const double x,
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


int gsl_sf_hyperg_1F1_impl(const double a, const double b, const double x,
                           double * result
                           )
{
  int a_neg_integer;    /*  a   negative integer  */
  int b_neg_integer;    /*  b   negative integer  */
  int bma_neg_integer;  /*  b-a negative integer  */

  double bma = b - a;
  double prec;

  a_neg_integer = ( a < 0.0  &&  fabs(a - rint(a)) < locEPS );
  b_neg_integer = ( b < 0.0  &&  fabs(b - rint(b)) < locEPS );
  bma_neg_integer = ( bma < 0.0  &&  fabs(bma - rint(bma)) < locEPS );

  /* case: a==b,  exp(x) */
  if(fabs(b-a) < locEPS) {
    return gsl_sf_exp_impl(x, result);
  }

  /* case: denominator zeroes before numerator */
  if(b_neg_integer && !(a_neg_integer && a > b + 0.1)) {
    return GSL_EDOM;
  }

  /* If a is a negative integer, then the
   * series truncates to a polynomial.
   */
  if(a_neg_integer) {
    double prec;
    return hyperg_1F1_series(a, b, x, result, &prec);
  }

  /* If b-a is a negative integer, use the Kummer transformation
   *    1F1(a,b,x) = Exp(x) 1F1(b-a,b,x)
   * to reduce it to a polynomial times an exponential.
   * Note that there can be no error condition here, since
   * there can only be an error if 'b' is a negative integer, but
   * in that case we would not have gotten this far unless 'a' was
   * a negative integer as well, in which case the above block
   * handled the situation.
   */
  if(bma_neg_integer) {
    double Ex, Kummer_1F1;
    int stat_E = gsl_sf_exp_impl(x, &Ex);
    int stat_K = hyperg_1F1_series(bma, b, -x, &Kummer_1F1, &prec);
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
  }

  
  /* Now we have dealt with any special negative integer cases,
   * including the error cases, so we are left with a well-behaved
   * series evaluation, though the arguments may be large.
   */
   
  
  
  
  
}
