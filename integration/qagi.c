#include <config.h>
#include <stdlib.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include <gsl_integration.h>

/* Evaluate an integral over an infinite range using the transformation

   integrate(f(x),-Inf,Inf) = integrate((f((1-t)/t) + f(-(1-t)/t))/t^2,0,1)

   */

static double transform (double t, void *params);

int
gsl_integration_qagi (gsl_function * f,
		      double epsabs, double epsrel,
		      gsl_integration_workspace * workspace,
		      size_t * last,
		      double *result, double *abserr, size_t * neval)
{
  int status;
  size_t nqeval = 0;
  gsl_integration_rule_t *integration_rule = &gsl_integration_qk15;

  gsl_function f_transform;

  f_transform.function = &transform;
  f_transform.params = f;

  status = gsl_integration_qags_impl (&f_transform, 0.0, 1.0, epsabs, epsrel,
				      workspace,
				      result, abserr, last, &nqeval,
				      integration_rule);

  /* convert from quadrature rule evaluations to function evaluations */

  *neval = 2 * 15 * nqeval;

  return status;
}

static double 
transform (double t, void *params)
{
  gsl_function *f = (gsl_function *) params;
  double x = (1 - t) / t;
  double y = GSL_FN_EVAL (f, x) + GSL_FN_EVAL (f, -x);
  double weight = t * t;
  return y / weight;
}
