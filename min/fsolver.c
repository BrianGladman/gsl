#include <config.h>
#include <stdlib.h>
#include <string.h>
#include <gsl_errno.h>
#include <gsl_roots.h>

gsl_min_fsolver *
gsl_min_fsolver_alloc (const gsl_min_fsolver_type * T, 
			 gsl_function * f, double minimum, gsl_interval x)
{
  int status;

  gsl_min_fsolver * s = (gsl_min_fsolver *) malloc (sizeof (gsl_min_fsolver));

  if (s == 0)
    {
      GSL_ERROR_RETURN ("failed to allocate space for min solver struct",
			GSL_ENOMEM, 0);
    };

  s->state = malloc (T->size);

  if (s->state == 0)
    {
      free (s);		/* exception in constructor, avoid memory leak */

      GSL_ERROR_RETURN ("failed to allocate space for min solver state",
			GSL_ENOMEM, 0);
    };

  s->type = T ;

  status = gsl_min_fsolver_set (s, f, minimum, x); /* seed the generator */

  if (status != GSL_SUCCESS)
    {
      free (s->state);
      free (s);		/* exception in constructor, avoid memory leak */

      GSL_ERROR_RETURN ("failed to set solver", status, 0);
    };

  return s;
}

int
gsl_min_fsolver_set (gsl_min_fsolver * s, gsl_function * f, double minimum, gsl_interval x)
{
  s->function = f;
  s->minimum = minimum;
  s->interval = x;

  if (x.lower > x.upper)
    {
      GSL_ERROR ("invalid interval (lower > upper)", GSL_EINVAL);
    }

  if (minimum > x.upper || minimum < x.lower) 
    {
      GSL_ERROR ("minimum must lie inside interval", GSL_EINVAL);
    }

  return (s->type->set) (s->state, s->function, &(s->minimum), &(s->interval));
}

int
gsl_min_fsolver_iterate (gsl_min_fsolver * s)
{
  return (s->type->iterate) (s->state, 
			     s->function, &(s->minimum), &(s->interval));
}

void
gsl_min_fsolver_free (gsl_min_fsolver * s)
{
  free (s->state);
  free (s);
}

const char *
gsl_min_fsolver_name (const gsl_min_fsolver * s)
{
  return s->type->name;
}

double
gsl_min_fsolver_minimum (const gsl_min_fsolver * s)
{
  return s->minimum;
}

gsl_interval
gsl_min_fsolver_interval (const gsl_min_fsolver * s)
{
  return s->interval;
}

