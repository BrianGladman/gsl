#include <config.h>
#include <stdlib.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include <gsl_integration.h>

/* Evaluate an integral over an infinite range using the transformation

   integrate(f(x),-Inf,b) = integrate(f(b-(1-t)/t)/t^2,0,1)

   */

struct params { double b ; gsl_function * f ; } ;

static double transform (double t, void *params);

int
gsl_integration_qagil (gsl_function * f,
		       double b,
		       double epsabs, double epsrel,
		       gsl_integration_workspace * workspace,
		       size_t * last,
		       double *result, double *abserr, size_t * neval)
{
  int status;
  size_t nqeval = 0;
  gsl_integration_rule_t *integration_rule = &gsl_integration_qk15;

  gsl_function f_transform;
  struct params transform_params  ;

  transform_params.b = b ;
  transform_params.f = f ;

  f_transform.function = &transform;
  f_transform.params = &transform_params;

  status = gsl_integration_qags_impl (&f_transform, 0.0, 1.0, epsabs, epsrel,
				      workspace,
				      result, abserr, last, &nqeval,
				      integration_rule);

  /* convert from quadrature rule evaluations to function evaluations */

  *neval = 15 * nqeval;

  return status;
}

static double 
transform (double t, void *params)
{
  struct params *p = (struct params *) params;
  double b = p->b;
  gsl_function * f = p->f;
  double x = b - (1 - t) / t;
  double y = GSL_FN_EVAL (f, x);
  double weight = t * t;
  return y / weight;
}
