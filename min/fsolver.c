#include <config.h>
#include <stdlib.h>
#include <string.h>
#include <gsl_errno.h>
#include <gsl_min.h>

gsl_min_fminimizer *
gsl_min_fminimizer_alloc (const gsl_min_fminimizer_type * T, 
			 gsl_function * f, double minimum, gsl_interval x)
{
  int status;

  gsl_min_fminimizer * s = (gsl_min_fminimizer *) malloc (sizeof (gsl_min_fminimizer));

  if (s == 0)
    {
      GSL_ERROR_RETURN ("failed to allocate space for minimizer struct",
			GSL_ENOMEM, 0);
    };

  s->state = malloc (T->size);

  if (s->state == 0)
    {
      free (s);		/* exception in constructor, avoid memory leak */

      GSL_ERROR_RETURN ("failed to allocate space for minimizer state",
			GSL_ENOMEM, 0);
    };

  s->type = T ;

  status = gsl_min_fminimizer_set (s, f, minimum, x); /* seed the generator */

  if (status != GSL_SUCCESS)
    {
      free (s->state);
      free (s);		/* exception in constructor, avoid memory leak */

      GSL_ERROR_RETURN ("failed to set minimizer", status, 0);
    };

  return s;
}

int
gsl_min_fminimizer_set (gsl_min_fminimizer * s, gsl_function * f, double minimum, gsl_interval x)
{
  s->function = f;
  s->minimum = minimum;
  s->interval = x;

  if (x.lower > x.upper)
    {
      GSL_ERROR ("invalid interval (lower > upper)", GSL_EINVAL);
    }

  if (minimum >= x.upper || minimum <= x.lower) 
    {
      GSL_ERROR ("minimum must lie inside interval (lower < x < upper)",
                 GSL_EINVAL);
    }

  return (s->type->set) (s->state, s->function, &(s->minimum), &(s->interval));
}

int
gsl_min_fminimizer_iterate (gsl_min_fminimizer * s)
{
  return (s->type->iterate) (s->state, 
			     s->function, &(s->minimum), &(s->interval));
}

void
gsl_min_fminimizer_free (gsl_min_fminimizer * s)
{
  free (s->state);
  free (s);
}

const char *
gsl_min_fminimizer_name (const gsl_min_fminimizer * s)
{
  return s->type->name;
}

double
gsl_min_fminimizer_minimum (const gsl_min_fminimizer * s)
{
  return s->minimum;
}

gsl_interval
gsl_min_fminimizer_interval (const gsl_min_fminimizer * s)
{
  return s->interval;
}

