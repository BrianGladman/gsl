/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include "gsl_sf_laguerre.h"
#include <gsl_errno.h>


double gsl_sf_laguerre_1(double a, double x)
{
  return 1. + a - x;
}

double gsl_sf_laguerre_2(double a, double x)
{
  double c0 = 0.5 * (2.+a)*(1.+a);
  double c1 = -2.-a;
  double c2 = -0.5/(2.+a);
  return c0 + c1*x*(1. + c2*x);
}

double gsl_sf_laguerre_3(double a, double x)
{
  double c0 = (3.+a)*(2.+a)*(1.+a) / 6.;
  double c1 = -c0 * 3. / (1.+a);
  double c2 = -1./(2.+a);
  double c3 = -1./(3.*(3.+a));
  return c0 + c1*x*(1. + c2*x*(1. + c3*x));
}

double gsl_sf_laguerre_4(double a, double x)
{
  double c0 = (4.+a)*(3.+a)*(2.+a)*(1.+a) / 24.;
  double c1 = -c0 * 4. / (1.+a);
  double c2 = -3./(2.*(2.+a));
  double c3 = -2./(3.*(3.+a));
  double c4 = -1./(4.*(4.+a));
  return c0 + c1*x*(1. + c2*x*(1. + c3*x*(1. + c4*x)));
}

double gsl_sf_laguerre_5(double a, double x)
{
  double c0 = (5.+a)*(4.+a)*(3.+a)*(2.+a)*(1.+a) / 120.;
  double c1 = -c0 * 5. / (1.+a);
  double c2 = -2./(2.+a);
  double c3 = -1./(3.+a);
  double c4 = -1./(2.*(4.+a));
  double c5 = -1./(5.*(5.+a));
  return c0 + c1*x*(1. + c2*x*(1. + c3*x*(1. + c4*x*(1. + c5*x))));
}

double gsl_sf_laguerre_cp(int n, double a, double x)
{
  if(n < 0) {
    char buff[100];
    sprintf(buff,"laguerre_cp: n=%d  < 0 not allowed", n);
    GSL_MESSAGE(buff);
    return 1.;
  }
  else if(n == 0) {
    return 1.;
  }
  else if(n == 1) {
    return 1. + a - x;
  }
  else {
    int i;
    double result;
    double * c = (double *)malloc((n+1) * sizeof(double));
    if(c == 0) {
      GSL_MESSAGE("gsl_sf_laguerre_cp: out of memory");
      return 0.;
    }
    else {
      c[0] = 1.;
      for(i=n; i>=1; i--) { c[0] *= (i + a)/i; }
  
      c[1] = -c[0] * n / (1.+a);
  
      for(i=2; i<=n; i++) { c[i] = -(n+1.-i)/(i*(i+a)); }

      result = 1. + c[n] * x;
      for(i=n-1; i>0; i--) {
        result = 1. + c[i] * x * result;
      }
      result += c[0] - 1.;

      free(c);
      return result;
    }
  }
}
