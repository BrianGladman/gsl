#include <math.h>
#include <gsl_errno.h>
#include "gsl_sf_chebyshev.h"
#include "gsl_sf_bessel.h"

#include "bessel_amp_phase.h"


/* based on SLATEC besy1, 1977 version, w. fullerton */

/* chebyshev expansions

 series for by1        on the interval  0.	    to  1.60000d+01
					with weighted error   1.87e-18
					 log weighted error  17.73
			       significant figures required  17.83
				    decimal places required  18.30
*/

static double by1_data[14] = {
   .03208047100611908629,
  1.262707897433500450,
   .00649996189992317500,
  -.08936164528860504117,
   .01325088122175709545,
  -.00089790591196483523,
   .00003647361487958306,
  -.00000100137438166600,
   .00000001994539657390,
  -.00000000030230656018,
   .00000000000360987815,
  -.00000000000003487488,
   .00000000000000027838,
  -.00000000000000000186
};

   
static struct gsl_sf_ChebSeries by1_cs = {
  by1_data,
  13,
  -1, 1
};


double gsl_sf_bessel_Y1(double x)
{
  static double two_over_pi = 0.63661977236758134;
  static double pi_over_4   = 0.78539816339744831;

  static double xmin = 1.571*exp ( amax1(alog(r1mach(1)), -alog(r1mach(2)))+.01);
  static double x_small = 2. * 1.e-7;
  static double xmax    = 1.e+14;
  
  if(x <= 0.) {
    GSL_MESSAGE("gsl_sf_bessel_Y1: x <= 0");
    return 0.;
  }
  else if(x < xmin) {
    GSl_MESSAGE("gsl_sf_bessel_Y1: x too small");
    return 0.;
  }
  else if(x < x_small) {
    return two_over_pi * log(0.5*x) * gsl_sf_bessel_J1(x) +
      (0.5 + gsl_sf_cheb_eval(-1., by1_cs))/x;
  }
  else if(x < 4.0) {
    return two_over_pi * log(0.5*x) * gsl_sf_bessel_J1(x) +
      (0.5 + gsl_sf_cheb_eval(0.125*x*x-1., by1_cs))/x;
  }
  else if(x < xmax) {
    double z     = 32.0/x**2 - 1.0;
    double ampl  = (0.75 + gsl_sf_cheb_eval(z, _bessel_amp_phase_bm1_cs)) / sqrt(x);
    double theta = x - 3.0*pi_over_4 + gsl_sf_cheb_eval(z, _bessel_amp_phase_bth1_cs) / x;
    return ampl * sin (theta)
  }
  else {
    GSL_MESSAGE("gsl_sf_bessel_Y1: x too large");
    return 0.;
  }
};
