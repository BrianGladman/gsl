#include <config.h>
#include <stdlib.h>
#include <math.h>
#include <gsl_math.h>
#include <gsl_integration.h>

struct fn_cauchy_params
{
  gsl_function *function;
  double singularity;
};


static double fn_cauchy (double t, void *params);

static void compute_moments (double cc,
			     double *cheb12, double *cheb24,
			     double *result, double *abserr);

void
qc25c (gsl_function * f, double a, double b, double c, 
       double *result, double *abserr);

void
qc25c (gsl_function * f, double a, double b, double c, 
       double *result, double *abserr)
{
  double cc = (2 * c - b - a) / (b - a);

  if (fabs (cc) > 1.1)
    {
      double resabs, resasc;

      gsl_function weighted_function;
      struct fn_cauchy_params fn_params;

      fn_params.function = f;
      fn_params.singularity = c;

      weighted_function.function = &fn_cauchy;
      weighted_function.params = &fn_params;

      gsl_integration_qk15 (&weighted_function, a, b, result, abserr,
			    &resabs, &resasc);
      return;
    }
  else
    {
      double cheb12[13], cheb24[25];
      gsl_integration_qcheb (f, a, b, cheb12, cheb24);
      compute_moments (cc, cheb12, cheb24, result, abserr);
      return;
    }
}

static double
fn_cauchy (double x, void *params)
{
  struct fn_cauchy_params *p = (struct fn_cauchy_params *) params;
  gsl_function *f = p->function;
  double c = p->singularity;
  return GSL_FN_EVAL (f, x) / (x - c);
}

static void
compute_moments (double cc, double *cheb12, double *cheb24,
		 double *result, double *abserr)
{
  size_t k;

  double amom0 = log (fabs ((1.0 - cc) / (1.0 + cc)));
  double amom1 = 2 + amom0 * cc;

  double res12 = cheb12[0] * amom0 + cheb12[1] * amom1;
  double res24 = cheb24[0] * amom0 + cheb24[1] * amom1;

  for (k = 2; k < 13; k++)
    {
      double amom2 = 2.0 * cc * amom1 - amom0;

      if ((k % 2) == 0)
	{
	  amom2 -= 4.0 / ((k - 1.0) * (k - 1.0) - 1.0);
	}

      res12 += cheb12[k] * amom2;
      res24 += cheb24[k] * amom2;

      amom0 = amom1;
      amom1 = amom2;
    }

  for (k = 13; k < 25; k++)
    {
      double amom2 = 2.0 * cc * amom1 - amom0;

      if ((k % 2) == 0)
	{
	  amom2 -= 4.0 / ((k - 1.0) * (k - 1.0) - 1.0);
	}

      res24 += cheb24[k] * amom2;

      amom0 = amom1;
      amom1 = amom2;
    }

  *result = res24;
  *abserr = fabs (res24 - res12);
}
