#include <math.h>
#include <gsl_errno.h>
#include "gsl_sf_cheb.h"
#include "gsl_sf_bessel.h"


/* based on SLATEC bk0(), bk0e() */

/* chebyshev expansions 

 series for bk0        on the interval  0.	    to  4.00000d+00
					with weighted error   3.57e-19
					 log weighted error  18.45
			       significant figures required  17.99
				    decimal places required  18.97

 series for ak0        on the interval  1.25000d-01 to  5.00000d-01
					with weighted error   5.34e-17
					 log weighted error  16.27
			       significant figures required  14.92
				    decimal places required  16.89

 series for ak02       on the interval  0.	    to  1.25000d-01
					with weighted error   2.34e-17
					 log weighted error  16.63
			       significant figures required  14.67
				    decimal places required  17.20
*/

static double bk0_data[11] = {
  -.03532739323390276872,
   .3442898999246284869, 
   .03597993651536150163,
   .00126461541144692592,
   .00002286212103119451,
   .00000025347910790261,
   .00000000190451637722,
   .00000000001034969525,
   .00000000000004259816,
   .00000000000000013744,
   .00000000000000000035
};

static struct gsl_sf_ChebSeries bk0_cs = {
  bk0_data,
  10,
  -1, 1
};

static double ak0_data[17] = {
  -.07643947903327941,
  -.02235652605699819,
   .00077341811546938,
  -.00004281006688886,
   .00000308170017386,
  -.00000026393672220,
   .00000002563713036,
  -.00000000274270554,
   .00000000031694296,
  -.00000000003902353,
   .00000000000506804,
  -.00000000000068895,
   .00000000000009744,
  -.00000000000001427,
   .00000000000000215,
  -.00000000000000033,
   .00000000000000005
};

static struct gsl_sf_ChebSeries ak0_cs = {
  ak0_data,
  16,
  -1, 1
};

static double ak02_data[14] = {
  -.01201869826307592,
  -.00917485269102569,
   .00014445509317750,
  -.00000401361417543,
   .00000015678318108,
  -.00000000777011043,
   .00000000046111825,
  -.00000000003158592,
   .00000000000243501,
  -.00000000000020743,
   .00000000000001925,
  -.00000000000000192,
   .00000000000000020,
  -.00000000000000002
};

static struct gsl_sf_ChebSeries ak02_cs = {
  ak02_data,
  13,
  -1, 1
};


double gsl_sf_bessel_K0_scaled(double x)
{
  static double x_small = 2.0 * 1.e-7;

  if(x <= 0.0) {
    GSL_MESSAGE("gsl_sf_bessel_K0_scaled: x <= 0");
    return 0.;
  }
  else if(x < x_small) {
    return  exp(x) * (-log(0.5*x)*gsl_sf_bessel_I0(x) - .25
      	      	      + gsl_sf_cheb_eval(-1., bk0_cs)
      	      	     )
  }
  else if(x <= 2.) {
    double y = x*x;
    return  exp(x) * (-log(0.5*x)*gsl_sf_bessel_I0(x) - .25
      	      	      + gsl_sf_cheb_eval(0.5*y-1., bk0_cs)
      	      	     )
  }
  else if(x <= 8.) {
    return (1.25 + gsl_sf_cheb_eval((16./x-5.)/3., ak0_cs)) / sqrt(x);
  }
  else {
    return (1.25 + gsl_sf_cheb_eval(16./x-1., ak02_cs)) / sqrt(x);
  } 
}


double gsl_sf_bessel_K0(double x)
{
     data ntk0, xsml, xmax / 0, 0., 0. /
c
  static double x_small = 2.*1.e-7;
  static double xmax ;
    /*
      xmax = -alog(r1mach(1))
      xmax = xmax - 0.5*xmax*alog(xmax)/(xmax+0.5) - 0.01
    */

  if(x <= 0.) {
    GSL_MESSAGE("gsl_sf_bessel_K0: x <= 0");
    return 0.;
  }
  else if(x < x_small) {
    return -log(0.5*x)*gsl_sf_bessel_I0(x) - .25 + gsl_sf_cheb_eval(-1., bk0_cs);
  }
  else if(x <= 2.) {
    double y = x*x;
    return -log(0.5*x)*gsl_sf_bessel_I0(x) - .25 + gsl_sf_cheb_eval(0.5*y-1., bk0_cs);
  }
  else if(x < xmax) {
    return exp(-x) * gsl_sf_bessel_K0_scaled(x);
  }
  else {
    return 0.; /* underflow ?? */
  }
}
