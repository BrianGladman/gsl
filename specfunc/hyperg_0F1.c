/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_exp.h"
#include "gsl_sf_gamma.h"
#include "gsl_sf_hyperg.h"

#define locEPS        (1000.0*GSL_MACH_EPS)
#define locMAX(a,b)   ((a) > (b) ? (a) : (b))
#define locMIN(a,b)   ((a) < (b) ? (a) : (b))


static
int
hyperg_0F1_series(double c, double x, double * result, double * prec)
{
  double cn  = c;
  double n   = 1.0;
  double sum = 1.0;
  double del = 1.0;
  double abs_del = 1.0;
  double max_abs_del = 1.0;
  double err;
  int n_min = 20;
  int n_max = 200;
  int n_far = 0;
  double special_term = 0.0;


  /* Figure out if there is a large contribution
   * hiding far out in the sum because c is
   * near a negative integer.
   * The constant to which we compare tells how
   * many significant figures can be lost when
   * near termination if we do not do some
   * special handling.
   */
  if(c < 0.0 && fabs(c - rint(c)) < 0.001) {
    double ln_term;
    double ln_poch;
    double ln_fact;
    int stat_exp;
    int stat_poch;
    double sign = (x < 0.0 && GSL_IS_ODD(n_far) ? -1.0 : 1.0);
    n_far = rint(c);
    gsl_sf_lnfact_impl(n_far, &ln_fact);
    stat_poch = gsl_sf_lnpoch_impl(c, n_far, &ln_poch);
    if(stat_poch != GSL_SUCCESS) {
      /* pochammer probably blew up, so there is no big contribution */
      special_term = 0.0;
      n_far = 0;
    }
    else {
      ln_term   = n_far * log(fabs(x)) - ln_fact - ln_poch;
      stat_exp  = gsl_sf_exp_impl(ln_term, &special_term);
      if(stat_exp == GSL_SUCCESS) {
        special_term *= sign;
      }
      else if(stat_exp == GSL_EUNDRFLW) {
        special_term = 0.0;
        n_far = 0;
      }
      else {
        *result = 0.0;
        return stat_exp;
      }
    }
  }

  while(abs_del/fabs(sum) > GSL_MACH_EPS || n < n_min) {
    double u, abs_u;

    if(cn == 0.0) {
      *result = 0.0;
      return GSL_EDOM;
    }
    if(n > n_max) {
      max_abs_del *= GSL_MACH_EPS;
      err     = fabs(GSL_MACH_EPS * n + max_abs_del);
      *prec   = err/(err + fabs(sum));
      *result = sum;
      if(*prec > locEPS)
        return GSL_ELOSS;
      else
        return GSL_SUCCESS;
    }
    
    u = x/(cn*n);
    abs_u = fabs(u);
    if(abs_u > 1.0 && max_abs_del > DBL_MAX/abs_u) {
      *prec   = 1.0;
      *result = sum;
      return GSL_ELOSS;
    }
    del *= u;

    sum += del;
    abs_del = fabs(del);
    max_abs_del = locMAX(abs_del, max_abs_del);

    cn += 1.0;
    n  += 1.0;
  }

  /* If the sum stopped before getting to the
   * distant large contribution, then include it now.
   */
  if(n_far > rint(n)) {
    sum += special_term;
  }

  max_abs_del *= GSL_MACH_EPS;
  err     = fabs(GSL_MACH_EPS * n + max_abs_del);
  *prec   = err/(err + fabs(sum));
  *result = sum;
  if(*prec > locEPS)
    return GSL_ELOSS;
  else
    return GSL_SUCCESS;
}



int
test_hyperg0F1_stuff(void)
{
  double c = -2.0+1.e-6;
  double x = 1.0;
  double result;
  double prec;
  int stat = hyperg_0F1_series(c, x, &result, &prec);
  
  printf("%20.16g  %10.6g    %20.16g  %20.16g  %s",
         c, x, result, prec, gsl_strerror(stat)
	 );
  printf("\n");
}


int
gsl_sf_hyperg_0F1_impl(double c, double x, double * result)
{
  int c_neg_integer = (c < 0.0 && fabs(c - rint(c)) < locEPS);

  if(fabs(c) < locEPS || c_neg_integer) {
    *result = 0.0;
    return GSL_EDOM;
  }

  if(x < 0.0) {
    double Jcm1;
    double lg_c;
    int stat_J = gsl_sf_bessel_Jnu();
  }
  else if(x == 0.0) {
    *result = 1.0;
    return GSL_SUCCESS;
  }
  else {
  }
}
