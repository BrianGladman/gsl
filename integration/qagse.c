#include <stdlib.h>
#include <gsl_errno.h>
#include <gsl_integration.h>

int
gsl_integration_qagse (double (*f)(double x),
		      double a, double b,
		      double epsabs, double epsrel,
		      size_t limit,
		      double alist[], double blist[], double rlist[], 
		      double elist[], size_t iord[], size_t * last,
		      double * result, double * abserr, size_t * neval)
{
  int status ;
  size_t nqeval = 0;
  gsl_integration_rule_t * integration_rule = &gsl_integration_qk21 ;

  status = gsl_integration_qagse_impl (f, a, b, epsabs, epsrel, limit,
				       result, abserr,
				       alist, blist, rlist, elist, iord, last,
				       &nqeval, 
				       integration_rule) ;

  /* convert from number of quadrature rule evaluations to number of
     function evaluations */

  *neval = 21 * nqeval ;

  return status ;
}

