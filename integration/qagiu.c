#include <config.h>
#include <stdlib.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include <gsl_integration.h>

#include "integration.h"

/* Evaluate an integral over an infinite range using the transformation

   integrate(f(x),a,Inf) = integrate(f(a+(1-t)/t)/t^2,0,1)

   */

struct params { double a ; gsl_function * f ; } ;

static double transform (double t, void *params);

int
gsl_integration_qagiu (gsl_function * f,
		       double a,
		       double epsabs, double epsrel, size_t limit,
		       gsl_integration_workspace * workspace,
		       double *result, double *abserr)
{
  int status;

  gsl_function f_transform;
  struct params transform_params  ;

  transform_params.a = a ;
  transform_params.f = f ;

  f_transform.function = &transform;
  f_transform.params = &transform_params;

  status = gsl_integration_qags_impl (&f_transform, 0.0, 1.0, 
				      epsabs, epsrel, limit,
				      workspace,
				      result, abserr,
				      &gsl_integration_qk15);

  return status;
}

static double 
transform (double t, void *params)
{
  struct params *p = (struct params *) params;
  double a = p->a;
  gsl_function * f = p->f;
  double x = a + (1 - t) / t;
  double y = GSL_FN_EVAL (f, x);
  return (y / t) / t;
}
