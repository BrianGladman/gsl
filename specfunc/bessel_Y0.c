#include <math.h>
#include <gsl_math.h>
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
  -.011277839392865573,
  -.12834523756042035,
  -.10437884799794249,
   .023662749183969695,
  -.002090391647700486,
   .000103975453939057,
  -.000003369747162423,
   .000000077293842676,
  -.000000001324976772,
   .000000000017648232,
  -.000000000000188105,
   .000000000000001641,
  -.000000000000000011
};

static struct gsl_sf_ChebSeries by0_cs = {
  by0_data,
  12,
  -1, 1
};


double gsl_sf_bessel_Y0(double x)
{
  const double two_over_pi = 2./M_PI;
  const double ln_half     = -M_LN2;
  const double x_small     = 2. * GSL_SQRT_MACH_EPS;
  const double xmax        = 1./GSL_MACH_EPS;

  if (x <= 0.) {
    GSL_ERROR_MESSAGE("gsl_sf_bessel_Y0: x <= 0", GSL_EDOM);
    return 0.;
  }
  else if(x < x_small){
    return two_over_pi*(ln_half + log(x))*gsl_sf_bessel_J0(x)
	    + .375 + gsl_sf_cheb_eval(-1., &by0_cs);
  }
  else if(x < 4.0) {
    return two_over_pi*(ln_half + log(x))*gsl_sf_bessel_J0(x)
	    + .375 + gsl_sf_cheb_eval(.125*x*x-1., &by0_cs);
  }
  else if(x < xmax) {
    double z     = 32.0/(x*x) - 1.0;
    double ampl  = (0.75 + gsl_sf_cheb_eval(z, &_bessel_amp_phase_bm0_cs)) / sqrt(x);
    double theta = x - M_PI_4 + gsl_sf_cheb_eval(z, &_bessel_amp_phase_bth0_cs) / x;
    return ampl * sin (theta);
  }
  else {
    GSL_ERROR_MESSAGE("gsl_sf_bessel_Y0: x too large", GSL_EUNDRFLW);
    return 0.;
  }
};
