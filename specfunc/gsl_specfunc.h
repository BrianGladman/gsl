/* $Id$ */
#ifndef _GSL_SPECFUNC_H_
#define _GSL_SPECFUNC_H_


/* Integer powers */
#define gsl_specfunc_SQR(x) ((x)*(x))
#define gsl_specfunc_CUBE(x) ((x)*(x)*(x))
double gsl_specfunc_pow_int(double x, int n);


/* Complex logarithm  (returns arg in [-pi,pi]) */
void gsl_specfunc_complex_log(double zr, double zi, double * lnr, double * theta);


/* Error function and associated functions. */
double gsl_specfunc_erf(double x);
double gsl_specfunc_erfc(double x);
double gsl_specfunc_erfc_asymptotic(double x);
double gsl_specfunc_erfseries(double x);
double gsl_specfunc_log_erfc(double x);
double gsl_specfunc_log_erfc_asymptotic(double x);
double gsl_specfunc_Q(double x);
double gsl_specfunc_Z(double x);


/* Chebyshev functions. */
struct gsl_specfunc_ChebSeries {
  double * c;   /* coefficients */
  int order;    /* order of expansion */
  double a;     /* lower interval point */
  double b;     /* upper interval point */
};
struct gsl_specfunc_ChebSeries * cheb_calc(double (*func)(double),
					double a, double b,
					int order);
double gsl_specfunc_cheb_eval(double x, const struct ChebSeries * cs);
void   gsl_specfunc_cheb_free(struct ChebSeries * cs);


/* Airy functions. */
double gsl_specfunc_Ai(double x);
double gsl_specfunc_Bi(double x);
double gsl_specfunc_Bi_scaled(double x); /* exponential prefactor removed when x > 0 */


#endif /* _GSL_SPECFUNC_H_ */
