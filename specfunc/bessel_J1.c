/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_chebyshev.h"
#include "gsl_sf_bessel.h"

#include "bessel_amp_phase.h"


/* based on SLATEC besj1, 1983 version, w. fullerton */

/* chebyshev expansions

 series for bj1        on the interval  0.	    to  1.60000d+01
					with weighted error   4.48e-17
					 log weighted error  16.35
			       significant figures required  15.77
				    decimal places required  16.89

*/

static double bj1_data[12] = {
  -.11726141513332787,
  -.25361521830790640,
   .050127080984469569,
  -.004631514809625081,
   .000247996229415914,
  -.000008678948686278,
   .000000214293917143,
  -.000000003936093079,
   .000000000055911823,
  -.000000000000632761,
   .000000000000005840,
  -.000000000000000044,
};
      

static struct gsl_sf_ChebSeries bj1_cs = {
  bj1_data,
  11,
  -1, 1
};

  
#define ROOT_EIGHT 2.82842712474619

double gsl_sf_bessel_J1(double x)
{ 
  const double x_small = ROOT_EIGHT * GSL_SQRT_MACH_EPS;
  const double xmax = 1./GSL_MACH_EPS;
  const double xmin = DBL_MIN;
  int err_status = GSL_SUCCESS;

  double y = abs(x);
  if(y < xmin) {
    err_status = GSL_EUNDRFLW;
    return 0.;
  }
  else if(y < x_small) {
    return 0.5*x;
  }
  else if(y < 4.0) {
    return x * (.25 + gsl_sf_cheb_eval(.125*y*y-1., &bj1_cs));
  }
  else if(y < xmax) {
    double z     = 32.0/(y*y) - 1.0;
    double ampl  = (0.75 + gsl_sf_cheb_eval(z, &_bessel_amp_phase_bm1_cs)) / sqrt(y);
    double theta = y - 3.0*M_PI_4 + gsl_sf_cheb_eval(z, &_bessel_amp_phase_bth1_cs) / y;
    return fortran_sign(ampl, x) * cos (theta);
  }
  else {
    GSL_ERROR_RETURN("gsl_sf_bessel_J1: x too large", GSL_EOVRFLW, 0.);
  }
}
