/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_gamma.h"
#include "gsl_sf_pow_int.h"
#include "gsl_sf_zeta.h"
#include "gsl_sf_fermi_dirac.h"

#define locEPS  (1000.0*GSL_MACH_EPS)


#ifndef Sqr
#define Sqr(x) ((x)*(x))
#endif


/* Handle case of integer j <= -2.
 */
static
int
fd_nint(const int j, const double x, double * result)
{
  const int nmax = 100;
  double qcoeff[nmax+1];

  if(j >= -1) {
    *result = 0.0;
    return GSL_ESANITY;
  }
  else if(j < -(nmax+1)) {
    *result = 0.0;
    return GSL_EUNIMPL;
  }
  else {
    double a, p, f;
    int i, k;
    int n = -(j+1);

    qcoeff[1] = 1.0;

    for(k=2; k<=n; k++) {
      qcoeff[k] = -qcoeff[k-1];
      for(i=k-1; i>=2; i--) {
        qcoeff[i] = i*qcoeff[i] - (k-(i-1))*qcoeff[i-1];
      }
    }

    if(x >= 0.0) {
      a = exp(-x);
      f = qcoeff[1];
      for(i=2; i<=n; i++) {
        f = f*a + qcoeff[i];
      }
    }
    else {
      a = exp(x);
      f = qcoeff[n];
      for(i=n-1; i>=1; i--) {
        f = f*a + qcoeff[i];
      }
    }

    p = gsl_sf_pow_int(1.0+a, j);
    *result = f*a*p;
  }
}


/* x < 0
 */
static
int
fd_neg(const double j, const double x, double * result)
{
  const int itmax = 100;
  double qnum[itmax+1], qden[itmax+1];

  if(x < GSL_LOG_DBL_MIN) {
    *result = 0.0;
    return GSL_SUCCESS;
  }
  else {
    double s;
    double xn = x;
    double ex  = -exp(x);
    double enx = -ex;
    double f = 0.0;
    int jterm;
    for(jterm=1; jterm<=itmax; jterm++) {
      double p = gsl_sf_pow_int(jterm, j+1);
      double f_previous = f;
      double term = enx/p;
      int stat_lu;
      xn += x;
      if(fabs(f-f_previous) < fabs(f)*10.0*GSL_MACH_EPS || xn < GSL_LOG_DBL_MIN) break;
      stat_lu = gsl_sf_levin_u_transform(term, jterm, qnum, qden, &f, s);
      enx *= ex;
    }
    if(jterm == itmax)
      return GSL_EMAXITER;
    else
      return GSL_SUCCESS;
  }
}


/* asymptotic expansion
 * j + 2.0 > 0.0
 */
static
int
fd_asymp(const double j, const double x, double * result)
{
  const int j_integer = ( fabs(j - floor(j+0.5)) < 100.0*GSL_MACH_EPS );
  const int itmax = 100;
  double lg;
  int stat_lg = gsl_sf_lngamma_impl(j + 2.0, &lg);
  double seqn = 0.5;
  double xm2  = (1.0/x)/x;
  double xgam = 1.0;
  double add  = DBL_MAX;
  double fneg;
  int stat_fneg;
  int stat_ser;
  int stat_eta = 0;
  int n;
  for(n=1; n<=itmax; n++) {
    double add_previous = add;
    double eta;
    gsl_sf_eta_int_impl(2*n, &eta);
    xgam = xgam * xm2 * (j + 1.0 - (2*n-2)) * (j + 1.0 - (2*n-1));
    add  = eta * xgam;
    if(fabs(add) > fabs(add_previous) && !j_integer) break;
    seqn += add;
  }
  stat_ser = ( fabs(add) > locEPS*fabs(seqn) ? GSL_ELOSS : GSL_SUCCESS );

  stat_fneg = fd_neg(j, -x, &fneg);
  *result = cos(j*M_PI) * fneg + 2.0 * seqn * exp((j+1.0)*log(x) - lg);
  return GSL_ERROR_SELECT_3(stat_fneg, stat_ser, stat_lg);
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


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

/* [Goano, TOMS-745, (4)] */
int gsl_sf_fermi_dirac_m1_impl(const double x, double * result)
{
  if(x < GSL_LOG_DBL_MIN) {
    *result = 0.0;
    return GSL_EUNDRFLW;
  }
  else if( x < 0.0) {
    double ex = exp(x);
    *result = ex/(1.0+ex);
    return GSL_SUCCESS;
  }
  else {
    *result = 1.0/(1.0 + exp(-x));
    return GSL_SUCCESS;
  }
}


/* [Goano, TOMS-745, (3)] */
int gsl_sf_fermi_dirac_0_impl(const double x, double * result)
{
  if(x < GSL_LOG_DBL_MIN) {
    *result = 0.0;
    return GSL_EUNDRFLW;
  }
  else if(x < -10) {
    double ex = exp(x);
    *result = ex * ( 1.0 - 0.5*ex + ex*ex/3.0 - ex*ex*ex/4.0);
    return GSL_SUCCESS;
  }
  else if(x < 10) {
    *result = log(1.0 + exp(x));
    return GSL_SUCCESS;
  }
  else {
    double ex = exp(-x);
    *result = x + ex * (1.0 - 0.5*ex + ex*ex/3.0 - ex*ex*ex/4.0);
    return GSL_SUCCESS;
  }
}


int gsl_sf_fermi_dirac_int_impl(const int j, const double x, double * result)
{
  if(j == 0) {
    return gsl_sf_fermi_dirac_0_impl(x, result);
  }
  else if(j == -1) {
    return gsl_sf_fermi_dirac_m1_impl(x, result);
  }
  else if(j < 0) {
    return fd_nint(j, x, result);
  }
  else if(x <= 0.0) {
    return fd_neg(j, x, result);
  }
  else {
    double k_div = -log10(10.0*GSL_MACH_EPS);
    double a1    = 2.0*k_div - j*(2.0+k_div/10.0);
    double a2    = sqrt(fabs((2.0*k_div - 1.0 - j)*(2.0*k_div - j)));
    double a     = locMAX(a1, a2);
    double xasymp = locMAX(j-1.0, a);
    double fasymp;
    int stat_asymp = fd_asymp(j, x, &fasymp);
  }
}


/* [Goano, TOMS-745, p. 222] */
int gsl_sf_fermi_dirac_inc_0_impl(const double x, const double b, double * result)
{
  if(b < 0.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else {
    double arg = b - x;
    double f0;
    int status = gsl_sf_fermi_integral_0_impl(arg, &f0);
    *result = f0 - arg;
    return status;
  }
}




int gsl_sf_fermi_integral_1_impl(const double A, double * result)
{
  if(A <= 0.) {
    return GSL_EDOM;
  }
  else if(A < 0.01) {
    /* asymptotic */
    double LA = -log(A);
    double rLA = sqrt(LA);
    double term1 = 2./3. * rLA * LA;
    double term2 = Sqr(M_PI)/(12. * rLA);
    double term3 = 7.*Sqr(M_PI)*Sqr(M_PI)/(960. * (rLA * LA*LA));
    *result = term1 + term2 + term3;
    return GSL_SUCCESS;
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
    *result = term1 + i1 + i2;
    return GSL_SUCCESS;
  }
  else if(A < 1.01) {
    /* Taylor series near A=1 */
    *result = 0.678148 * (1 - 0.791 * (A-1.));
    return GSL_SUCCESS;
  }
  else if(A < 20.) {
    /* integrate */
    double integral =
      trap_rule_2(F1_c_integrand, 0., 4., A, 50)
      + trap_rule_2(F1_c_integrand, 4., 6., A, 50)
      + trap_rule_2(F1_c_integrand, 6., 10., A, 50)
      + trap_rule_2(F1_c_integrand, 10., 20., A, 50);
    *result = integral / A;
    return GSL_SUCCESS;
  }
  else {
    /* asymptotic */
    double coeff = 0.5 * sqrt(M_PI);
    double sum = 0.;
    int i;
    for(i=1; i<4; i++) { sum += 1./(i*sqrt(i)) / gsl_sf_pow_int(-A, i); }
    *result = coeff / A * (1. + sum);
    return GSL_SUCCESS;
  }
}

int gsl_sf_fermi_integral_2_impl(const double A, double * result)
{
  if(A <= 0.) {
    return GSL_EDOM;
  }
  else if(A < 0.01) {
    /* asymptotic */
    double LA = -log(A);
    double rLA = sqrt(LA);
    double term1 = rLA;
    double term2 = -Sqr(M_PI)/24. / (rLA * LA);
    double term3 = -7. * Sqr(M_PI)*Sqr(M_PI) / 384. / (rLA * LA*LA*LA);
    *result = term1 + term2 + term3;
    return GSL_SUCCESS;
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

    *result = term1 + term2;
    return GSL_SUCCESS;
  }
  else if(A < 1.01) {
    /* Taylor series near A=1 */
    *result = 0.536077 * (1 - 0.466 * (A-1.));
    return GSL_SUCCESS;
  }
  else if(A < 20.) {
    /* integrate */
    double integral =
      trap_rule_2(F2_c_integrand, 0., 4., A, 100)
      + trap_rule_2(F2_c_integrand, 4., 6., A, 100)
      + trap_rule_2(F2_c_integrand, 6., 10., A, 50)
      + trap_rule_2(F2_c_integrand, 10., 20., A, 50);
    *result = integral / A;
    return GSL_SUCCESS;
  }
  else {
    /* asymptotic */
    double coeff = 0.5 * sqrt(M_PI);
    double sum = 0.;
    int i;
    for(i=1; i<4; i++) { sum += sqrt(1./(double)i) / gsl_sf_pow_int(-A, i); }
    *result = coeff / A * (1. + sum);
    return GSL_SUCCESS;
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_fermi_dirac_0_e(const double x, double * result)
{
  int status = gsl_sf_fermi_dirac_0_impl(x, result);
  if(status != GSL_SUCCESS){
    GSL_ERROR("gsl_sf_fermi_dirac_0_e", status);
  }
  return status;
}

int gsl_sf_fermi_dirac_m1_e(const double x, double * result)
{
  int status = gsl_sf_fermi_dirac_m1_impl(x, result);
  if(status != GSL_SUCCESS){
    GSL_ERROR("gsl_sf_fermi_dirac_m1_e", status);
  }
  return status;
}


int 
gsl_sf_fermi_dirac_inc_0_e(const double x, const double b, double * result)
{
  int status = gsl_sf_fermi_dirac_inc_0_impl(x, b, result);
  if(status != GSL_SUCCESS){
    GSL_ERROR("gsl_sf_fermi_dirac_inc_0_e", status);
  }
  return status;
}


int gsl_sf_fermi_integral_1_e(const double A, double * result)
{
  int status = gsl_sf_fermi_integral_1_impl(A, result);
  if(status != GSL_SUCCESS){
    GSL_ERROR("gsl_sf_fermi_integral_1_e", status);
  }
  return status;
}

int gsl_sf_fermi_integral_2_e(const double A, double * result)
{
  int status = gsl_sf_fermi_integral_2_impl(A, result);
  if(status != GSL_SUCCESS){
    GSL_ERROR("gsl_sf_fermi_integral_2_e", status);
  }
  return status;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*-*/

double gsl_sf_fermi_dirac_0(const double x)
{
  double y;
  int status = gsl_sf_fermi_dirac_0_impl(x, &y);
  if(status != GSL_SUCCESS){
    GSL_WARNING("gsl_sf_fermi_dirac_0", status);
  }
  return y;
}

double gsl_sf_fermi_dirac_m1(const double x)
{
  double y;
  int status = gsl_sf_fermi_dirac_m1_impl(x, &y);
  if(status != GSL_SUCCESS){
    GSL_WARNING("gsl_sf_fermi_dirac_m1", status);
  }
  return y;
}


double gsl_sf_fermi_dirac_inc_0(const double x, const double b)
{
  double y;
  int status = gsl_sf_fermi_dirac_inc_0_impl(x, b, &y);
  if(status != GSL_SUCCESS){
    GSL_WARNING("gsl_sf_fermi_dirac_inc_0", status);
  }
  return y;
}


double gsl_sf_fermi_integral_1(const double A)
{
  double y;
  int status = gsl_sf_fermi_integral_1_impl(A, &y);
  if(status != GSL_SUCCESS){
    GSL_WARNING("gsl_sf_fermi_integral_1", status);
  }
  return y;
}

double gsl_sf_fermi_integral_2(const double A)
{
  double y;
  int status = gsl_sf_fermi_integral_2_impl(A, &y);
  if(status != GSL_SUCCESS){
    GSL_WARNING("gsl_sf_fermi_integral_2", status);
  }
  return y;
}

double gsl_sf_fermi_zeta_inverse(const double dbar, const double g, const double prec)
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

double gsl_sf_fermi_dirac(const double beta, const double zeta_inverse, const double E)
{
  return 1./(zeta_inverse*exp(beta*E)+1.);
}
