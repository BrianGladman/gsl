/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_chebyshev.h"
#include "gsl_sf_gamma.h"
#include "gsl_sf_pow_int.h"
#include "gsl_sf_zeta.h"
#include "gsl_sf_fermi_dirac.h"


#define locMAX(a,b)   ((a) > (b) ? (a) : (b))
#define locEPS        (1000.0*GSL_MACH_EPS)


/* Chebyshev fit for F_{1/2}(3/4(t+1) + 1/2);  -1 < t < 1, 1/2 < x < 2
 */
static double fd_half_data[20] = {
  3.817142729374783,
  0.85335973724616,
  0.06209005531895366,
 -0.0001736543402433676,
 -0.0001353392641073591,
  8.89368871503926e-6,
  2.986477928734743e-7,
 -8.20194437034516e-8,
  3.025685657576104e-9,
  4.701755984548583e-10,
 -5.773975741973913e-11,
 -2.016720124231596e-13,
  5.413336445769799e-13,
 -3.491651412446117e-14,
 -3.69704267200177e-15,
  1.443289932012703e-16,
  2.775557561562891e-17,
  1.84297022087776e-15,
 -1.076916333886401e-15,
  2.636779683484746e-16
};
static struct gsl_sf_cheb_series fd_half_cs = {
  fd_half_data,
  19,
  -1, 1,
  (double *)0,
  (double *)0
};



/* Goano's modification of the Levin-u implementation from TOMS-602.
 */
static
int
fd_whiz(const double term, const int iterm,
        double * qnum, double * qden,
        double * result, double * s)
{
  if(iterm == 0) *s = 0.0;

  *s += term;

  qden[iterm] = 1.0/(term*(iterm+1.0)*(iterm+1.0));
  qnum[iterm] = *s * qden[iterm];

  if(iterm > 0) {
    double factor = 1.0;
    double ratio  = iterm/(iterm+1.0);
    int j;
    for(j=iterm-1; j>=0; j--) {
      double c = factor * (j+1.0) / (iterm+1.0);
      factor *= ratio;
      qden[j] = qden[j+1] - c * qden[j];
      qnum[j] = qnum[j+1] - c * qnum[j];
    }
  }

  *result = qnum[0] / qden[0];
  return GSL_SUCCESS;
}


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
    return GSL_SUCCESS;
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
    for(jterm=0; jterm<=itmax; jterm++) {
      double p = gsl_sf_pow_int(jterm+1, j+1);
      double f_previous = f;
      double term = enx/p;
      fd_whiz(term, jterm, qnum, qden, &f, &s);
      xn += x;
      if(fabs(f-f_previous) < fabs(f)*10.0*GSL_MACH_EPS || xn < GSL_LOG_DBL_MIN) break;
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
  const int itmax = 200;
  double lg;
  int stat_lg = gsl_sf_lngamma_impl(j + 2.0, &lg);
  double seqn = 0.5;
  double xm2  = (1.0/x)/x;
  double xgam = 1.0;
  double add = DBL_MAX;
  double fneg;
  int stat_fneg;
  int stat_ser;
  int n;
  for(n=1; n<=itmax; n++) {
    double add_previous = add;
    double eta;
    gsl_sf_eta_int_impl(2*n, &eta);
    xgam = xgam * xm2 * (j + 1.0 - (2*n-2)) * (j + 1.0 - (2*n-1));
    add  = eta * xgam;
    if(!j_integer && fabs(add) > fabs(add_previous)) break;
    if( j_integer && 2*n > j) break;
    seqn += add;
  }
  stat_ser = ( fabs(add) > locEPS*fabs(seqn) ? GSL_ELOSS : GSL_SUCCESS );

  stat_fneg = fd_neg(j, -x, &fneg);
  *result = cos(j*M_PI) * fneg + 2.0 * seqn * exp((j+1.0)*log(x) - lg);
  return GSL_ERROR_SELECT_3(stat_ser, stat_fneg, stat_lg);
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
  else if(x < -5.0) {
    double ex  = exp(x);
    double ser = 1.0 - ex*(0.5 - ex*(1.0/3.0 - ex*(1.0/4.0 - ex*(1.0/5.0 - ex/6.0))));
    *result = ex * ser;
    return GSL_SUCCESS;
  }
  else if(x < 10.0) {
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


int gsl_sf_fermi_dirac_half_impl(const double x, double * result)
{
  if(x < GSL_LOG_DBL_MIN) {
    *result = 0.0;
    return GSL_EUNDRFLW;
  }
  else if(x < -0.5) {
    /* series [Goano (6)]
     */
    double ex   = exp(x);
    double term = ex;
    double sum  = term;
    int n;
    for(n=2; n<100 ; n++) {
      term *= -ex * pow((n-1.0)/n, 1.5);
      sum  += term;
      if(fabs(term/sum) < GSL_MACH_EPS) break;
    }
    *result = sum;
    return GSL_SUCCESS;
  }
  else if(x < 0.5) {
    /* economized Taylor series
     */
    const double n1 = 6.218242902223359 + x;
    const double n2 = 16.541445090841950 - 0.2793612233580265*x + x*x;
    const double n3 = 41.317548401295120 + 10.448057133628715*x + x*x;
    const double d1 = 58.802032649535580 - 13.415748084443560*x + x*x;
    const double d2 = 23.780193903718020 - 3.0871842455825022*x + x*x;
    const double d3 = 12.538103971158665 - 0.4487609507636724*x + x*x;
    const double num = n1*n2*n3;
    const double den = d1*d2*d3;
    *result = 3.15652196749738490 * num / den;
    return GSL_SUCCESS;
  }
  else if(x < 2.0) {
    /* Chebyshev fit.
     */
    double t = 4.0/3.0*(x - 0.5) - 1.0;
    *result  = gsl_sf_cheb_eval(&fd_half_cs, t);
    return GSL_SUCCESS;
  }
  else {
    return fd_asymp(0.5, x, result);
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
    int status = gsl_sf_fermi_dirac_0_impl(arg, &f0);
    *result = f0 - arg;
    return status;
  }
}



/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_fermi_dirac_m1_e(const double x, double * result)
{
  int status = gsl_sf_fermi_dirac_m1_impl(x, result);
  if(status != GSL_SUCCESS){
    GSL_ERROR("gsl_sf_fermi_dirac_m1_e", status);
  }
  return status;
}


int gsl_sf_fermi_dirac_0_e(const double x, double * result)
{
  int status = gsl_sf_fermi_dirac_0_impl(x, result);
  if(status != GSL_SUCCESS){
    GSL_ERROR("gsl_sf_fermi_dirac_0_e", status);
  }
  return status;
}


int gsl_sf_fermi_dirac_int_e(const int j, const double x, double * result)
{
  int status = gsl_sf_fermi_dirac_int_impl(j, x, result);
  if(status != GSL_SUCCESS){
    GSL_ERROR("gsl_sf_fermi_dirac_int_e", status);
  }
  return status;
}


int gsl_sf_fermi_dirac_half_e(const double x, double * result)
{
  int status = gsl_sf_fermi_dirac_half_impl(x, result);
  if(status != GSL_SUCCESS){
    GSL_ERROR("gsl_sf_fermi_dirac_half_e", status);
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


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*-*/

double gsl_sf_fermi_dirac_m1(const double x)
{
  double y;
  int status = gsl_sf_fermi_dirac_m1_impl(x, &y);
  if(status != GSL_SUCCESS){
    GSL_WARNING("gsl_sf_fermi_dirac_m1", status);
  }
  return y;
}


double gsl_sf_fermi_dirac_0(const double x)
{
  double y;
  int status = gsl_sf_fermi_dirac_0_impl(x, &y);
  if(status != GSL_SUCCESS){
    GSL_WARNING("gsl_sf_fermi_dirac_0", status);
  }
  return y;
}


double gsl_sf_fermi_dirac_int(const int j, const double x)
{
  double y;
  int status = gsl_sf_fermi_dirac_int_impl(j, x, &y);
  if(status != GSL_SUCCESS){
    GSL_WARNING("gsl_sf_fermi_dirac_int", status);
  }
  return y;
}


double gsl_sf_fermi_dirac_half(const double x)
{
  double y;
  int status = gsl_sf_fermi_dirac_half_impl(x, &y);
  if(status != GSL_SUCCESS){
    GSL_WARNING("gsl_sf_fermi_dirac_half", status);
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
