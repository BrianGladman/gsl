/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_exp.h"
#include "gsl_sf_hyperg.h"


int
gsl_sf_hyperg_0F1_impl(double c, double x, double * result)
{
  if(x < 0.0) {
    double Jcm1;
    double lg_c;
    int stat_J = gsl_sf_bessel_Jnu();
  }
  else if(x == 0.0) {
    *result = 1.0;
    return GSL_SUCCESS;
  }
  else {
  }
}

