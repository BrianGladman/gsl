#include <config.h>
#include <stdlib.h>
#include <gsl_errno.h>
#include <gsl_integration.h>

#include "integration.h"

int
gsl_integration_qags (const gsl_function *f,
		      double a, double b,
		      double epsabs, double epsrel, size_t limit,
		      gsl_integration_workspace * workspace,
		      double * result, double * abserr)
{
  int status = gsl_integration_qags_impl (f, a, b, epsabs, epsrel, limit,
					  workspace, 
					  result, abserr, 
					  &gsl_integration_qk21) ;
  return status ;
}

