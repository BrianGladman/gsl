#include <config.h>
#include <stdlib.h>
#include <string.h>
#include <gsl_errno.h>
#include <gsl_roots.h>

gsl_root_fdf_solver *
gsl_root_fdf_solver_alloc (const gsl_root_fdf_solver_type * T, 
			   gsl_fdf * f, double root)
{

  gsl_root_fdf_solver * s = (gsl_root_fdf_solver *) malloc (sizeof (gsl_root_fdf_solver));

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

  gsl_root_fdf_solver_set (s, f, root); /* seed the generator */

  return s;
}

int
gsl_root_fdf_solver_set (gsl_root_fdf_solver * s, gsl_fdf * f, double root)
{
  s->fdf = f;
  s->root = root;

  return (s->set) (s->state, s->fdf, &(s->root));
}

int
gsl_root_fdf_solver_iterate (gsl_root_fdf_solver * s)
{
  return (s->iterate) (s->state, s->fdf, &(s->root));
}

void
gsl_root_fdf_solver_free (gsl_root_fdf_solver * s)
{
  free (s->state);
  free (s);
}

