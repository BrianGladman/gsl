#include <config.h>
#include <stdlib.h>
#include <string.h>
#include <gsl_errno.h>
#include <gsl_roots.h>

gsl_root_fdfsolver *
gsl_root_fdfsolver_alloc (const gsl_root_fdfsolver_type * T, 
			   gsl_function_fdf * f, double root)
{

  gsl_root_fdfsolver * s = (gsl_root_fdfsolver *) malloc (sizeof (gsl_root_fdfsolver));

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

  s->type = T ;

  gsl_root_fdfsolver_set (s, f, root); /* seed the generator */

  return s;
}

int
gsl_root_fdfsolver_set (gsl_root_fdfsolver * s, gsl_function_fdf * f, double root)
{
  s->fdf = f;
  s->root = root;

  return (s->type->set) (s->state, s->fdf, &(s->root));
}

int
gsl_root_fdfsolver_iterate (gsl_root_fdfsolver * s)
{
  return (s->type->iterate) (s->state, s->fdf, &(s->root));
}

void
gsl_root_fdfsolver_free (gsl_root_fdfsolver * s)
{
  free (s->state);
  free (s);
}

const char *
gsl_root_fdfsolver_name (const gsl_root_fdfsolver * s)
{
  return s->type->name;
}

double
gsl_root_fdfsolver_root (const gsl_root_fdfsolver * s)
{
  return s->root;
}


