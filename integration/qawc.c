#include <config.h>
#include <stdlib.h>
#include <gsl_errno.h>
#include <gsl_integration.h>

int
gsl_integration_qawc (const gsl_function *f,
		      double a, double b, double c,
		      double epsabs, double epsrel,
		      gsl_integration_workspace * workspace,
		      size_t * last,
		      double * result, double * abserr, size_t * neval)
{
  int status ;
  size_t nqeval = 0;

  status = gsl_integration_qawc_impl (f, a, b, c, epsabs, epsrel,
				      workspace, 
				      result, abserr, last, &nqeval) ;
  
  /* convert from quadrature rule evaluations to function evaluations */

  *neval = 21 * nqeval ; /* FIXME */

  return status ;
}

