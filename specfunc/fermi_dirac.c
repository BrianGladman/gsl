/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_fermi_dirac.h"

#ifndef Sqr
#define Sqr(x) ((x)*(x))
#endif


double gsl_sf_fermi_dirac(double beta, double zeta_inverse, double E)
{
  return 1./(zeta_inverse*exp(beta*E)+1.);
}

/* implement trapezoid rule for func with 2 args
   integrates with respect to the first arg
   p = parameters
   nsteps = number of steps
 */
static double trap_rule_2(double(*f)(double, double),
			  double a, double b,
			  double p,
			  int nsteps)
{
  double delta = fabs(b - a) / (double) nsteps;
  double ends = 0.5 * delta * (f(a, p) + f(b, p));
  double sum = 0.;
  int i;
  for(i = 1; i<nsteps; i++) sum += f(a + i * delta, p);
  return ends + delta * sum;
}

/* F1 first integrand for A < 1 */
static double F1_a_integrand(double t, double abs_lnA)
{
  double term1 = sqrt(abs_lnA + t);
  double term2 = sqrt(abs_lnA - t);
  double numer = term1 - term2;
  double denom = 1. + exp(t);
  return numer / denom;
}

/* F1 second integrand for A < 1 */
static double F1_b_integrand(double t, double abs_lnA)
{
  double numer = sqrt(t + abs_lnA);
  double denom = 1. + exp(t);
  return numer / denom;
}

/* F1 integrand for A > 1 */
static double F1_c_integrand(double y, double A)
{
  double denom = exp(y) + 1/A;
  double numer = sqrt(y);
  return numer / denom;
}


double gsl_sf_fermi_integral_1(double A)
{
  if(A <= 0.) {
    char buff[100];
    sprintf(buff, "fermi_integral_1: bad arg  %18.12g\n", A);
    GSL_ERROR_RETURN(buff, GSL_EDOM, 6600.);
  }
  else if(A < 0.01) {
    /* asymptotic */
    double LA = -log(A);
    double rLA = sqrt(LA);
    double term1 = 2./3. * rLA * LA;
    double term2 = Sqr(M_PI)/(12. * rLA);
    double term3 = 7.*Sqr(M_PI)*Sqr(M_PI)/(960. * (rLA * LA*LA));
    return term1 + term2 + term3;
  }
  else if(A < 0.99) {
    /* integrate */
    double LA = -log(A);
    double i1 = trap_rule_2(F1_a_integrand, 0., LA, LA, 100);
    double i2 = 
        trap_rule_2(F1_b_integrand, LA,    LA+2.,  LA, 100)
      + trap_rule_2(F1_b_integrand, LA+2., LA+4.,  LA, 50)
      + trap_rule_2(F1_b_integrand, LA+4., LA+6.,  LA, 50)
      + trap_rule_2(F1_b_integrand, LA+6., LA+10., LA, 50);
    double term1 = 2./3. * LA * sqrt(LA);
    return term1 + i1 + i2;
  }
  else if(A < 1.01) {
    /* Taylor series near A=1 */
    return 0.678148 * (1 - 0.791 * (A-1.));
  }
  else if(A < 20.) {
    /* integrate */
    double integral =
      trap_rule_2(F1_c_integrand, 0., 4., A, 50)
      + trap_rule_2(F1_c_integrand, 4., 6., A, 50)
      + trap_rule_2(F1_c_integrand, 6., 10., A, 50)
      + trap_rule_2(F1_c_integrand, 10., 20., A, 50);
    return integral / A;
  }
  else {
    /* asymptotic */
    double coeff = 0.5 * sqrt(M_PI);
    double sum = 0.;
    int i;
    for(i=1; i<4; i++) { sum += pow(1./(double)i, 1.5) / pow(-A, i); }
    return coeff / A * (1. + sum);
  }
}

/* F2 first integrand */
static double F2_a_integrand(double t, double abs_lnA)
{
  double lt = log(t);
  double numer = sqrt(1. + lt/abs_lnA);
  double denom = (t + 1.)*(t + 1.);
  return numer/denom;
}

/* F2 second integrand */
static double F2_b_integrand(double t, double abs_lnA)
{
  double et = exp(-t);
  double numer = et * sqrt((1.+t/abs_lnA)/(1.+1./abs_lnA));
  double denom = (1.+et)*(1.+et);
  return numer/denom;
}

/* F2 third integrand */
static double F2_c_integrand(double y, double A)
{
  double ey = exp(-y);
  double numer = sqrt(y) * ey;
  double denom = (1.+ey/A)*(1.+ey/A);
  return numer/denom;
}

double gsl_sf_fermi_integral_2(double A)
{
  if(A <= 0.) {
    char buff[100];
    sprintf(buff, "fermi_integral_2: bad arg  %18.12g\n", A);
    GSL_ERROR_RETURN(buff, GSL_EDOM, 6666.);
    return 6666.;
  }
  else if(A < 0.01) {
    /* asymptotic */
    double LA = -log(A);
    double rLA = sqrt(LA);
    double term1 = rLA;
    double term2 = -Sqr(M_PI)/24. / (rLA * LA);
    double term3 = -7. * Sqr(M_PI)*Sqr(M_PI) / 384. / (rLA * LA*LA*LA);
    return term1 + term2 + term3;
  }
  else if(A < 0.99) {
    /* integrate */
    double LA = -log(A);
    double i1 = trap_rule_2(F2_a_integrand, A, exp(1.), LA, 150);
    double i2 = 
      trap_rule_2(F2_b_integrand, 1., 4., LA, 100)
      + trap_rule_2(F2_b_integrand, 4., 6.,LA, 50)
      + trap_rule_2(F2_b_integrand, 6., 10., LA, 50)
      + trap_rule_2(F2_b_integrand, 10., 20., LA, 50);

    double term1 = sqrt(LA) * i1;
    double term2 = sqrt(1.+LA) * i2;

    return term1 + term2;
  }
  else if(A < 1.01) {
    /* Taylor series near A=1 */
    return 0.536077 * (1 - 0.466 * (A-1.));
  }
  else if(A < 20.) {
    /* integrate */
    double integral =
      trap_rule_2(F2_c_integrand, 0., 4., A, 100)
      + trap_rule_2(F2_c_integrand, 4., 6., A, 100)
      + trap_rule_2(F2_c_integrand, 6., 10., A, 50)
      + trap_rule_2(F2_c_integrand, 10., 20., A, 50);
    return integral / A;
  }
  else {
    /* asymptotic */
    double coeff = 0.5 * sqrt(M_PI);
    double sum = 0.;
    int i;
    for(i=1; i<4; i++) { sum += sqrt(1./(double)i) / pow(-A, i); }
    return coeff / A * (1. + sum);
  }
}


double gsl_sf_fermi_zeta_inverse(double dbar, double g, double prec)
{
  double coeff = g / (4. * M_PI * M_PI);
  double zeta_inverse_0 = 1.;
  double zeta_inverse_1 = 1.;
  double guess;

  /* find a value where dbar is too high */
  guess = coeff * fermi_integral_1(zeta_inverse_0);
  while(guess < dbar) {
    zeta_inverse_0 *= 0.5;
    guess = coeff * fermi_integral_1(zeta_inverse_0);
  }

  /* find a value where dbar is too low */
  guess = coeff * fermi_integral_1(zeta_inverse_1);
  while(guess > dbar) {
    zeta_inverse_1 *= 2.;
    guess = coeff * fermi_integral_1(zeta_inverse_1);
  }

  /* the two values of zeta_inverse now bracket
     the solution; iterate the bracketing
   */
  while(1) {
    double new_zeta_inverse = 0.5 * (zeta_inverse_0 + zeta_inverse_1);
    guess = coeff * fermi_integral_1(new_zeta_inverse);

    if(fabs((guess - dbar)/dbar) < prec) {
      return new_zeta_inverse;
    }
    else if(guess > dbar) {
      zeta_inverse_0 = new_zeta_inverse;
    }
    else if(guess < dbar) {
      zeta_inverse_1 = new_zeta_inverse;
    }
  }
}
