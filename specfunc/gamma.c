/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <stdio.h>
#include <math.h>
#include <gsl_errno.h>
#include <gsl_math.h>
#include "gsl_sf_log.h"
#include "gsl_sf_trig.h"
#include "gsl_sf_gamma.h"

extern double hypot(double, double);


#define LogRootTwoPi_ 0.918938533204673
#define Max(a,b) ((a) > (b) ? (a) : (b))


static double gamma_cof[6] = {76.18009172947146, -86.50532032941677,
			    	  24.01409824083091, -1.231739572450155,
			    	  0.1208650973866179e-2, -0.5395239384953e-5
	    	    	    	  };

double gsl_sf_lngamma(double xx)
{
  int j;
  double x = xx;
  double y = xx;
  double tmp = x + 5.5 - (x + 0.5) * log(x + 5.5);
  double ser = 1.000000000190015;
  
  for(j=0; j<=5; j++) ser += gamma_cof[j]/(++y);
  
  return -tmp + log(2.5066282746310005 * ser / x);
}


/* series which appears in stirling approximation for ln(gamma(z)), z = x + iy;
   see Carrier,Krook,Pearson p.188
   the coefficients are c_n = B_2n/((2n-1)(2n))
 */
static void ln_stirling_series(double x, double y, double * sx, double * sy)
{
  int i;
  double r2 = x*x + y*y;
  double z_inv_pow_x[18], z_inv_pow_y[18];
  static double c[10] = 
    {1., 1./12., -1./360., 1./1260., -1./1680., 1./1188.,
     -691./360360., 1./156., -3617./122400., 43867./244188.
    };

  z_inv_pow_x[0] = 1.;
  z_inv_pow_y[0] = 0.;
  z_inv_pow_x[1] =  x/r2;
  z_inv_pow_y[1] = -y/r2;
  for(i=2; i<18; i++) { /* create powers of 1/z */
    z_inv_pow_x[i] = z_inv_pow_x[1]*z_inv_pow_x[i-1] - z_inv_pow_y[1]*z_inv_pow_y[i-1];
    z_inv_pow_y[i] = z_inv_pow_x[1]*z_inv_pow_y[i-1] + z_inv_pow_y[1]*z_inv_pow_x[i-1];
  }

  *sx = 0.;
  *sy = 0.;
  for(i=1; i<10; i++) {
    *sx += c[i] * z_inv_pow_x[2*i-1];
    *sy += c[i] * z_inv_pow_y[2*i-1];
  }
}


void gsl_sf_complex_lngamma(double zr, double zi, double * lnr, double * arg)
{
  double stirl_cut = 50.;        /* use stirling above this point */
  double za = hypot(zr, zi);
  double x, y, r, theta;         /* transformed variables */
  double t_lnr, t_arg;           /* log(abs) and arg of intermediate result */
  double lr;                     /* log(r) */
  double ln_sr, s_theta;         /* Stirling series result */

  /* normalization factor from recursion step */
  double ln_renorm_r, ln_renorm_theta;

  if(za == 0.) {
    char buff[100];
    sprintf(buff,"gsl_sf_complex_lngamma: zr= %g  zi= %g", zr, zi);
    GSL_ERROR(buff, GSL_EDOM);
    *lnr = 0.;
    *arg = 0.;
    return;
  }

  /* transform to positive real part using reflection */
  if(zr <= 0.) {
    x = 1.-zr;
    y = -zi;
  }
  else {
    x = zr;
    y = zi;
  }

  /* make (x,y) large and keep track of the
     renormalization factor from the recursion relation
   */
  ln_renorm_r = 0.;
  ln_renorm_theta = 0.;
  if(za < stirl_cut) {
    int n;
    int num_recurse = Max((int)(stirl_cut - x), 0);
    for(n=0; n<num_recurse; n++) {
      double a,b;
      gsl_sf_complex_log(x, y, &a, &b);
      ln_renorm_r += a;
      ln_renorm_theta += b;
      x += 1.;
    }
  }

  /* calculate Gamma(x+iy) using Stirling */
  r = hypot(x, y);
  theta = atan2(y, x);
  lr = log(r);
  ln_stirling_series(x, y, &ln_sr, &s_theta);
  t_lnr = LogRootTwoPi_ - x - theta*y + (x-0.5)*lr + ln_sr;
  t_arg = theta*x + y * lr - y + s_theta;

  /* divide by the renormalization factor */
  t_lnr -= ln_renorm_r;
  t_arg -= ln_renorm_theta;

  /* undo the reflection */
  if(zr <= 0.0) {
    double prefactor_x, prefactor_y;
    double denom_x, denom_y, denom_abs2;

    gsl_sf_complex_sin(M_PI*zr, M_PI*zi, &denom_x, &denom_y);
    denom_abs2 = denom_x*denom_x + denom_y*denom_y;

    if(denom_abs2 < 1.e-60) {
      char buff[100];
      sprintf(buff,"complex_lngamma: z= (%g,%g) near a negative integer",
	      zr, zi);
      GSL_ERROR(buff, GSL_EDOM);
      *lnr = 0.;
      *arg = 0.;
    }
    prefactor_x =  M_PI*denom_x/denom_abs2;
    prefactor_y = -M_PI*denom_y/denom_abs2;
    *lnr = 0.5*log(prefactor_x*prefactor_x + prefactor_y*prefactor_y) - t_lnr;
    *arg = -t_arg + atan2(prefactor_y, prefactor_x);
  }
  else {
    *lnr = t_lnr;
    *arg = t_arg;
  }

  gsl_sf_angle_restrict_symm(arg);
}

