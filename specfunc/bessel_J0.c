/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_chebyshev.h"
#include "gsl_sf_bessel.h"

#include "bessel_amp_phase.h"


/* based on SLATEC besj0, 1977 version, w. fullerton */

/* chebyshev expansions for Bessel functions

 series for bj0        on the interval  0.	    to  1.60000d+01
					with weighted error   7.47e-18
					 log weighted error  17.13
			       significant figures required  16.98
				    decimal places required  17.68

*/

static double bj0_data[13] = {
   .100254161968939137, 
  -.665223007764405132, 
   .248983703498281314, 
  -.0332527231700357697,
   .0023114179304694015,
  -.0000991127741995080,
   .0000028916708643998,
  -.0000000612108586630,
   .0000000009838650793,
  -.0000000000124235515,
   .0000000000001265433,
  -.0000000000000010619,
   .0000000000000000074,
};
  	

static struct gsl_sf_ChebSeries bj0_cs = {
  bj0_data,
  12,
  -1, 1
};


double gsl_sf_bessel_J0(double x)
{
  const double x_small = 2. * GSL_SQRT_MACH_EPS;
  const double xmax = 1./GSL_MACH_EPS;

  double y = fabs(x);

  if(y < x_small) {
    return 1.;
  }
  else if(y <= 4.0) {
    return gsl_sf_cheb_eval(0.125*y*y - 1.0, &bj0_cs);
  }
  else if (y < xmax) {
    double z     = 32.0/(y*y) - 1.0;
    double ampl  = (0.75 + gsl_sf_cheb_eval(z, &_bessel_amp_phase_bm0_cs)) / sqrt(y);
    double theta = y - M_PI_4 + gsl_sf_cheb_eval(z, &_bessel_amp_phase_bth0_cs) / y;
    return ampl * cos(theta);
  }
  else {
    GSL_ERROR_RETURN("gsl_sf_bessel_J0: |x| too large", GSL_EOVRFLW, 0.);
  }
}
