/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_errno.h>
#include "gsl_sf_log.h"


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_complex_log_impl(double zr, double zi, double * lnr, double * theta)
{
  if(zr != 0.0 || zi != 0.0) {
    double r2 = zr*zr + zi*zi;
    *lnr = 0.5*log(r2);
    *theta = atan2(zi, zr);
    return GSL_SUCCESS;
  }
  else {
    return GSL_EDOM;
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* Error Handling Versions *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_complex_log_e(double zr, double zi, double * lnr, double * theta)
{
  int status = gsl_sf_complex_log_impl(zr, zi, lnr, theta);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_complex_log_e", status);
  }
  return status;
}
