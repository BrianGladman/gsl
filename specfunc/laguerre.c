/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <gsl_errno.h>
#include "gsl_sf_laguerre.h"


/*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*/

/* polynomial based on confluent hypergeometric representation */
static int laguerre_n_cp(const int n, const double a, const double x, double * result)
{
  double * c = (double *)malloc((n+1) * sizeof(double));
  if(c == 0) {
    return GSL_ENOMEM;
  }
  else {
    int i;
    double ans;

    c[0] = 1.;
    for(i=n; i>=1; i--) { c[0] *= (i + a)/i; }

    c[1] = -c[0] * n / (1.+a);

    for(i=2; i<=n; i++) { c[i] = -(n+1.-i)/(i*(i+a)); }

    ans = 1. + c[n] * x;
    for(i=n-1; i>0; i--) {
      ans = 1. + c[i] * x * ans;
    }
    ans += c[0] - 1.;

    free(c);
    *result = ans;
    return GSL_SUCCESS;
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_laguerre_n_impl(const int n, const double a, const double x, double * result)
{
  if(n < 0 || a <= -1. || x < 0.) {
    return GSL_EDOM;
  }
  else if(n == 0) {
    *result = 1.;
    return GSL_SUCCESS;
  }
  else if(n == 1) {
    *result = 1. + a - x;
    return GSL_SUCCESS;
  }
  else {
    int status = laguerre_n_cp(n, a, x, result);
    return status;
  }
}

/*-*-*-*-*-*-*-*-*-*-*-* Error Handling Versions *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_laguerre_n_e(int n, double a, double x, double * result)
{
  int status = gsl_sf_laguerre_n_impl(n, a, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_laguerre_n_e", status);
  }
  return status;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*-*/

double gsl_sf_laguerre_1(const double a, const double x)
{
  return 1. + a - x;
}

double gsl_sf_laguerre_2(const double a, const double x)
{
  double c0 = 0.5 * (2.+a)*(1.+a);
  double c1 = -(2.+a);
  double c2 = -0.5/(2.+a);
  return c0 + c1*x*(1. + c2*x);
}

double gsl_sf_laguerre_3(const double a, const double x)
{
  double c0 = (3.+a)*(2.+a)*(1.+a) / 6.;
  double c1 = -c0 * 3. / (1.+a);
  double c2 = -1./(2.+a);
  double c3 = -1./(3.*(3.+a));
  return c0 + c1*x*(1. + c2*x*(1. + c3*x));
}

double gsl_sf_laguerre_4(const double a, const double x)
{
  double c0 = (4.+a)*(3.+a)*(2.+a)*(1.+a) / 24.;
  double c1 = -c0 * 4. / (1.+a);
  double c2 = -3./(2.*(2.+a));
  double c3 = -2./(3.*(3.+a));
  double c4 = -1./(4.*(4.+a));
  return c0 + c1*x*(1. + c2*x*(1. + c3*x*(1. + c4*x)));
}

double gsl_sf_laguerre_5(const double a, const double x)
{
  double c0 = (5.+a)*(4.+a)*(3.+a)*(2.+a)*(1.+a) / 120.;
  double c1 = -c0 * 5. / (1.+a);
  double c2 = -2./(2.+a);
  double c3 = -1./(3.+a);
  double c4 = -1./(2.*(4.+a));
  double c5 = -1./(5.*(5.+a));
  return c0 + c1*x*(1. + c2*x*(1. + c3*x*(1. + c4*x*(1. + c5*x))));
}

double gsl_sf_laguerre_n(int n, double a, double x)
{
  double y;
  int status = gsl_sf_laguerre_n_impl(n, a, x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_laguerre_n_e", status);
  }
  return y;
}
