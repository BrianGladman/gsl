/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_errno.h>
#include "gsl_sf_log.h"


int gsl_sf_complex_log(double zr, double zi, double * lnr, double * theta)
{
  if(zr != 0.0 || zi != 0.0) {
    double r2 = zr*zr + zi*zi;
    *lnr = 0.5*log(r2);
    *theta = atan2(zi, zr);
    return GSL_SUCCESS;
  }
  else {
    GSL_ERROR("gsl_sf_complex_log: z=0.0", GSL_EDOM);
  }
}
