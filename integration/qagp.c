#include <config.h>
#include <stdlib.h>
#include <gsl_errno.h>
#include <gsl_integration.h>

#include "integration.h"

int
gsl_integration_qagp (const gsl_function *f,
		      double * pts, size_t npts,
		      double epsabs, double epsrel, size_t limit,
		      gsl_integration_workspace * workspace,
		      double * result, double * abserr)
{
  int status = gsl_integration_qagp_impl (f, pts, npts, 
					  epsabs, epsrel, limit,
					  workspace,
					  result, abserr,
					  &gsl_integration_qk21) ;
  
  return status ;
}

