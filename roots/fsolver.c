#include <config.h>
#include <stdlib.h>
#include <string.h>
#include <gsl_errno.h>
#include <gsl_roots.h>

gsl_root_f_solver *
gsl_root_f_solver_alloc (const gsl_root_f_solver_type * T, 
			 gsl_function * f, gsl_interval x)
{
  int status;

  gsl_root_f_solver * s = (gsl_root_f_solver *) malloc (sizeof (gsl_root_f_solver));

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

  status = gsl_root_f_solver_set (s, f, x); /* seed the generator */

  if (status != GSL_SUCCESS)
    {
      free (s->state);
      free (s);		/* exception in constructor, avoid memory leak */

      GSL_ERROR_RETURN ("failed to set solver", status, 0);
    };

  return s;
}

int
gsl_root_f_solver_set (gsl_root_f_solver * s, gsl_function * f, gsl_interval x)
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
gsl_root_f_solver_iterate (gsl_root_f_solver * s)
{
  return (s->iterate) (s->state, 
			      s->function, &(s->root), &(s->interval));
}

void
gsl_root_f_solver_free (gsl_root_f_solver * s)
{
  free (s->state);
  free (s);
}

