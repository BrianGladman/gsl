#include <config.h>
#include <stdlib.h>
#include <string.h>
#include <gsl_errno.h>
#include <gsl_roots.h>

gsl_root_fsolver *
gsl_root_fsolver_alloc (const gsl_root_fsolver_type * T, 
			 gsl_function * f, gsl_interval x)
{
  int status;

  gsl_root_fsolver * s = (gsl_root_fsolver *) malloc (sizeof (gsl_root_fsolver));

  if (s == 0)
    {
      GSL_ERROR_RETURN ("failed to allocate space for root solver struct",
			GSL_ENOMEM, 0);
    };

  s->state = malloc (T->size);

  if (s->state == 0)
    {
      free (s);		/* exception in constructor, avoid memory leak */

      GSL_ERROR_RETURN ("failed to allocate space for root solver state",
			GSL_ENOMEM, 0);
    };

  s->name = T->name ;
  s->size = T->size ;
  s->set = T->set ;
  s->iterate = T->iterate ;

  status = gsl_root_fsolver_set (s, f, x); /* seed the generator */

  if (status != GSL_SUCCESS)
    {
      free (s->state);
      free (s);		/* exception in constructor, avoid memory leak */

      GSL_ERROR_RETURN ("failed to set solver", status, 0);
    };

  return s;
}

int
gsl_root_fsolver_set (gsl_root_fsolver * s, gsl_function * f, gsl_interval x)
{
  s->function = f;
  s->root = 0;
  s->interval = x;

  if (x.lower > x.upper)
    {
      GSL_ERROR ("invalid interval (lower > upper)", GSL_EINVAL);
    }

  return (s->set) (s->state, s->function, &(s->root), &(s->interval));
}

int
gsl_root_fsolver_iterate (gsl_root_fsolver * s)
{
  return (s->iterate) (s->state, 
			      s->function, &(s->root), &(s->interval));
}

void
gsl_root_fsolver_free (gsl_root_fsolver * s)
{
  free (s->state);
  free (s);
}

const char *
gsl_root_fsolver_name (const gsl_root_fsolver * s)
{
  return s->name;
}

double
gsl_root_fsolver_root (const gsl_root_fsolver * s)
{
  return s->root;
}

gsl_interval
gsl_root_fsolver_interval (const gsl_root_fsolver * s)
{
  return s->interval;
}


