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

static void compute_moments (double cc, double *moment);

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
      double cheb12[13], cheb24[25], moment[25];
      double res12 = 0, res24 = 0;
      size_t i;
      gsl_integration_qcheb (f, a, b, cheb12, cheb24);
      compute_moments (cc, moment);

      for (i = 0; i < 13; i++)
	{
	  res12 += cheb12[i] * moment[i];
	}

      for (i = 0; i < 25; i++)
	{
	  res24 += cheb24[i] * moment[i];
	}

      *result = res24;
      *abserr = fabs(res24 - res12) ;

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
compute_moments (double cc, double *moment)
{
  size_t k;

  double a0 = log (fabs ((1.0 - cc) / (1.0 + cc)));
  double a1 = 2 + a0 * cc;

  moment[0] = a0;
  moment[1] = a1;

  for (k = 2; k < 25; k++)
    {
      double a2;

      if ((k % 2) == 0)
	{
	  const double km1 = k - 1.0;
	  a2 = 2.0 * cc * a1 - a0 - 4.0 / (km1 * km1 - 1.0);
	}
      else
	{
	  a2 = 2.0 * cc * a1 - a0;
	}

      moment[k] = a2;

      a0 = a1;
      a1 = a2;
    }
}
