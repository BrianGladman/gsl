#include <config.h>
#include <stdlib.h>
#include <string.h>
#include <gsl_errno.h>
#include <gsl_multiroots.h>

gsl_multiroot_fsolver *
gsl_multiroot_fsolver_alloc (const gsl_multiroot_fsolver_type * T, 
                             gsl_multiroot_function * f, gsl_vector * x)
{
  int status;

  gsl_multiroot_fsolver * s;

  const size_t n = f->n;

  if (x->size != n) 
    {
      GSL_ERROR_RETURN ("vector length not compatible with function", 
                        GSL_EBADLEN, 0);
    }

  s = (gsl_multiroot_fsolver *) malloc (sizeof (gsl_multiroot_fsolver));

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

  status = gsl_multiroot_fsolver_set (s, f, x); /* seed the generator */

  if (status != GSL_SUCCESS)
    {
      free (s->state);
      free (s);		/* exception in constructor, avoid memory leak */

      GSL_ERROR_RETURN ("failed to set solver", status, 0);
    };

  return s;
}

int
gsl_multiroot_fsolver_set (gsl_multiroot_fsolver * s, gsl_function * f, gsl_vector * x)
{
  s->function = f;
  s->root = x;

  return (s->type->set) (s->state, s->function, s->root);
}

int
gsl_multiroot_fsolver_iterate (gsl_multiroot_fsolver * s)
{
  return (s->type->iterate) (s->state, s->function, s->root);
}

void
gsl_multiroot_fsolver_free (gsl_multiroot_fsolver * s)
{
  free (s->state);
  free (s);
}

const char *
gsl_multiroot_fsolver_name (const gsl_multiroot_fsolver * s)
{
  return s->type->name;
}

double
gsl_multiroot_fsolver_root (const gsl_multiroot_fsolver * s)
{
  return s->root;
}
