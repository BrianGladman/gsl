#include <config.h>
#include <stdlib.h>
#include <gsl_errno.h>
#include <gsl_integration.h>

#include "integration.h"

int
gsl_integration_qawc (gsl_function *f,
		      double a, double b, double c,
		      double epsabs, double epsrel, size_t limit,
		      gsl_integration_workspace * workspace,
		      double * result, double * abserr)
{
  int status = gsl_integration_qawc_impl (f, a, b, c, epsabs, epsrel, limit,
					  workspace, result, abserr) ;
  
  return status ;
}



