#include <math.h>
#include <gsl_errno.h>
#include "gsl_sf_chebyshev.h"
#include "gsl_sf_bessel.h"

#include "bessel_amp_phase.h"


/* based on SLATEC besy0, 1980 version, w. fullerton */

/* chebyshev expansions

 series for by0        on the interval  0.	    to  1.60000d+01
					with weighted error   1.20e-17
					 log weighted error  16.92
			       significant figures required  16.15
				    decimal places required  17.48
*/

static double by0_data[13] = {
  -.0112778393 92865573e0
  -.1283452375 6042035e0
  -.1043788479 9794249e0
   .0236627491 83969695e0
  -.0020903916 47700486e0
   .0001039754 53939057e0
  -.0000033697 47162423e0
   .0000000772 93842676e0
  -.0000000013 24976772e0
   .0000000000 17648232e0
  -.0000000000 00188105e0
   .0000000000 00001641e0
  -.0000000000 00000011e0
};

static struct gsl_sf_ChebSeries by0_cs = {
  by0_data,
  12,
  -1, 1
};


double gsl_sf_bessel_Y0(double x)
{
  static double two_over_pi = 0.63661977236758134;
  static double pi_over_4   = 0.78539816339744831;
  static double ln_half     = -0.693147180559945309;
  static double x_small     = 2. * 1.e-7;
  static double xmax        = 1.e+14;

  if (x <= 0.) {
    GSL_MESSAGE("gsl_sf_bessel_Y0: x <= 0");
    return 0.;
  }
  else if(x < x_small){
    return two_over_pi*(ln_half + log(x))*gsl_sf_bessel_J0(x)
	    + .375 + gsl_sf_cheb_eval(-1., by0_cs);
  }
  else if(x < 4.0) {
    return two_over_pi*(ln_half + log(x))*gsl_sf_bessel_J0(x)
	    + .375 + gsl_sf_cheb_eval(.125*x*x-1., by0_cs);
  }
  else if(x < xmax) {
    double z     = 32.0/x**2 - 1.0;
    double ampl  = (0.75 + gsl_sf_cheb_eval(z, _bessel_amp_phase_bm0_cs)) / sqrt(x);
    double theta = x - pi_over_4 + gsl_sf_cheb_eval(z, _bessel_amp_phase_bth0_cs) / x;
    return ampl * sin (theta)
  }
  else {
    GSL_MESSAGE("gsl_sf_bessel_Y0: x too large");
    return 0.;
  }
};
