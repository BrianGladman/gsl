/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_SF_POW_INT_H_
#define GSL_SF_POW_INT_H_

/* Calculate small integer powers.
 * Does not check for overflow/underflow.
 */
double gsl_sf_pow_2(const double x);
double gsl_sf_pow_3(const double x);
double gsl_sf_pow_4(const double x);
double gsl_sf_pow_5(const double x);
double gsl_sf_pow_6(const double x);
double gsl_sf_pow_7(const double x);
double gsl_sf_pow_8(const double x);
double gsl_sf_pow_9(const double x);


/* Calculate x^n.
 * Does not check for overflow/underflow.
 */
double gsl_sf_pow_int(double x, int n);


#ifdef HAVE_INLINE
extern inline double gsl_sf_pow_2(const double x) { return x*x;   }
extern inline double gsl_sf_pow_3(const double x) { return x*x*x; }
extern inline double gsl_sf_pow_4(const double x) { double x2 = x*x;   return x2*x2;	}
extern inline double gsl_sf_pow_5(const double x) { double x2 = x*x;   return x2*x2*x;  }
extern inline double gsl_sf_pow_6(const double x) { double x2 = x*x;   return x2*x2*x2; }
extern inline double gsl_sf_pow_7(const double x) { double x3 = x*x*x; return x3*x3*x;  }
extern inline double gsl_sf_pow_8(const double x) { double x2 = x*x;   double x4 = x2*x2; return x4*x4; }
extern inline double gsl_sf_pow_9(const double x) { double x3 = x*x*x; return x3*x3*x3; }
#endif /* HAVE_INLINE */


#endif /* GSL_SF_POW_INT_H_ */
