/* brent.c -- brent minimum finding algorithm */

#include <config.h>

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include <gsl_math.h>
#include <gsl_errno.h>
#include <gsl_min.h>

#include "min.h"

typedef struct
  {
    double d, e, v, w;
    double f_v, f_w;
    double f_lower, f_upper, f_minimum;
  }
brent_state_t;

int brent_init (void *vstate, gsl_function * f, double *minimum, gsl_interval * x);
int brent_iterate (void *vstate, gsl_function * f, double *minimum, gsl_interval * x);

int
brent_init (void *vstate, gsl_function * f, double *minimum, gsl_interval * x)
{
  brent_state_t *state = (brent_state_t *) vstate;

  const double golden = 0.318966;	/* golden = (3 - sqrt(5))/2 */

  double x_lower = x->lower;
  double x_upper = x->upper;

  double v = x_lower + golden * (x_upper - x_lower);
  double w = v;

  double f_lower, f_upper, f_minimum, f_vw;

  state->v = v;
  state->w = w;

  state->d = 0;
  state->e = 0;

  SAFE_FUNC_CALL (f, x_lower, &f_lower);
  SAFE_FUNC_CALL (f, x_upper, &f_upper);
  SAFE_FUNC_CALL (f, *minimum, &f_minimum);
  SAFE_FUNC_CALL (f, v, &f_vw);

  state->f_lower = f_lower;
  state->f_upper = f_upper;
  state->f_minimum = f_minimum;
  state->f_v = f_vw;
  state->f_w = f_vw;

  if (f_minimum >= f_lower || f_minimum >= f_upper)
    {
      GSL_ERROR ("endpoints do not enclose a minimum", GSL_EINVAL);
    }

  return GSL_SUCCESS;

}

int
brent_iterate (void *vstate, gsl_function * f, double *minimum, gsl_interval * x)
{
  brent_state_t *state = (brent_state_t *) vstate;

  const double x_lower = x->lower;
  const double x_upper = x->upper;

  const double z = *minimum;
  double d = state->e;
  double e = state->d;
  double u, f_u;
  const double v = state->v;
  const double w = state->w;
  const double f_v = state->f_v;
  const double f_w = state->f_w;
  const double f_z = state->f_minimum;

  const double golden = 0.318966;	/* golden = (3 - sqrt(5))/2 */

  const double w_lower = (z - x_lower);
  const double w_upper = (x_upper - z);

  const double tolerance =  GSL_SQRT_DBL_EPSILON * fabs (z);

  double p = 0, q = 0, r = 0;

  const double midpoint = 0.5 * (x_lower + x_upper);

  if (fabs (e) > tolerance)
    {
      /* fit parabola */

      r = (z - w) * (f_z - f_v);
      q = (z - v) * (f_z - f_w);
      p = (z - v) * q - (z - w) * r;
      q = 2 * (q - r);

      if (q > 0)
	{
	  p = -p;
	}
      else
	{
	  q = -q;
	}

      r = e;
      e = d;
    }

  if (fabs (p) < fabs (0.5 * q * r) && p < q * w_lower && p < q * w_upper)
    {
      double t2 = 2 * tolerance ;

      d = p / q;
      u = z + d;

      if ((u - x_lower) < t2 || (x_upper - z) < t2)
	{
	  d = (z < midpoint) ? tolerance : -tolerance ;
	}
    }
  else
    {
      e = (z < midpoint) ? x_upper - z : -(z - x_lower) ;
      d = golden * e;
    }


  if (fabs (d) >= tolerance)
    {
      u = z + d;
    }
  else
    {
      u = z + ((d > 0) ? tolerance : -tolerance) ;
    }

  state->e = e;
  state->d = d;

  SAFE_FUNC_CALL(f, u, &f_u);

  if (f_u > f_z)
    {
      if (u < z)
	{
	  x->lower = u;
          state->f_lower = f_u;
          return GSL_SUCCESS;
	}
      else
	{
	  x->upper = u;
          state->f_upper = f_u;
          return GSL_SUCCESS;
	}
    }
  else if (f_u < f_z)
    {
      if (u < z)
	{
	  x->upper = z;
          state->f_upper = f_z;
	}
      else
	{
	  x->lower = z;
          state->f_lower = f_z;
	}

      state->v = w;
      state->f_v = f_w;
      state->w = z;
      state->f_w = f_z;
      *minimum = u;
      state->f_minimum = f_u;
      return GSL_SUCCESS;
    }
  else if (f_u <= f_w || w == z)
    {
      state->v = w;
      state->f_v = f_w;
      state->w = u;
      state->f_w = f_u;
      return GSL_SUCCESS;
    }
  else if (f_u <= f_v || v == z || v == w)
    {
      state->v = u;
      state->f_v = f_u;
      return GSL_SUCCESS;
    }
  else
    {
      return GSL_FAILURE;
    }
}


static const gsl_min_fminimizer_type brent_type =
{"brent",			/* name */
 sizeof (brent_state_t),
 &brent_init,
 &brent_iterate};

const gsl_min_fminimizer_type *gsl_min_fminimizer_brent = &brent_type;
