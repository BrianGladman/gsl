#include <config.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_multiroots.h>

gsl_multiroot_fsolver *
gsl_multiroot_fsolver_alloc (const gsl_multiroot_fsolver_type * T, 
                             gsl_multiroot_function * f, gsl_vector * x)
{
  int status;

  gsl_multiroot_fsolver * s;

  const size_t n = f->n ;

  if (x->size != n) 
    {
      GSL_ERROR_RETURN ("vector length not compatible with function", 
                        GSL_EBADLEN, 0);
    }

  s = (gsl_multiroot_fsolver *) malloc (sizeof (gsl_multiroot_fsolver));

  if (s == 0)
    {
      GSL_ERROR_RETURN ("failed to allocate space for multiroot solver struct",
			GSL_ENOMEM, 0);
    }

  s->x = gsl_vector_calloc (n);

  if (s->x == 0) 
    {
      free (s);
      GSL_ERROR_RETURN ("failed to allocate space for x", GSL_ENOMEM, 0);
    }

  s->f = gsl_vector_calloc (n);

  if (s->f == 0) 
    {
      gsl_vector_free (s->x);
      free (s);
      GSL_ERROR_RETURN ("failed to allocate space for f", GSL_ENOMEM, 0);
    }

  s->dx = gsl_vector_calloc (n);

  if (s->dx == 0) 
    {
      gsl_vector_free (s->x);
      gsl_vector_free (s->f);
      free (s);
      GSL_ERROR_RETURN ("failed to allocate space for dx", GSL_ENOMEM, 0);
    }

  s->state = malloc (T->size);

  if (s->state == 0)
    {
      gsl_vector_free (s->dx);
      gsl_vector_free (s->x);
      gsl_vector_free (s->f);
      free (s);		/* exception in constructor, avoid memory leak */
      
      GSL_ERROR_RETURN ("failed to allocate space for multiroot solver state",
			GSL_ENOMEM, 0);
    }

  s->type = T ;

  status = (s->type->alloc)(s->state, n);

  if (status != GSL_SUCCESS)
    {
      (s->type->free)(s->state);
      free (s->state);
      gsl_vector_free (s->dx);
      gsl_vector_free (s->x);
      gsl_vector_free (s->f);
      free (s);		/* exception in constructor, avoid memory leak */
      
      GSL_ERROR_RETURN ("failed to set solver", status, 0);
    }
  
  status = gsl_multiroot_fsolver_set (s, f, x); /* seed the generator */
  
  if (status != GSL_SUCCESS)
    {
      free (s->state);
      gsl_vector_free (s->dx);
      gsl_vector_free (s->x);
      gsl_vector_free (s->f);
      free (s);		/* exception in constructor, avoid memory leak */
      
      GSL_ERROR_RETURN ("failed to set solver", status, 0);
    }

  return s;
}

int
gsl_multiroot_fsolver_set (gsl_multiroot_fsolver * s, 
                           gsl_multiroot_function * f, 
                           gsl_vector * x)
{
  s->function = f;
  gsl_vector_cpy(s->x,x);
  
  return (s->type->set) (s->state, s->function, s->x, s->f, s->dx);
}

int
gsl_multiroot_fsolver_iterate (gsl_multiroot_fsolver * s)
{
  return (s->type->iterate) (s->state, s->function, s->x, s->f, s->dx);
}

void
gsl_multiroot_fsolver_free (gsl_multiroot_fsolver * s)
{
  (s->type->free) (s->state);
  free (s->state);
  gsl_vector_free (s->dx);
  gsl_vector_free (s->x);
  gsl_vector_free (s->f);
  free (s);
}

const char *
gsl_multiroot_fsolver_name (const gsl_multiroot_fsolver * s)
{
  return s->type->name;
}

gsl_vector *
gsl_multiroot_fsolver_root (const gsl_multiroot_fsolver * s)
{
  return s->x;
}
