#include <config.h>
#include <stdlib.h>
#include <stdio.h>
#include <gsl_errno.h>
#include <gsl_rng.h>

gsl_rng *
gsl_rng_alloc (const gsl_rng_type * T)
{

  gsl_rng *r = (gsl_rng *) malloc (sizeof (gsl_rng));

  if (r == 0)
    {
      GSL_ERROR_RETURN ("failed to allocate space for rng struct",
			GSL_ENOMEM, 0);
    };

  r->state = malloc (T->size);

  if (r->state == 0)
    {
      free (r);		/* exception in constructor, avoid memory leak */

      GSL_ERROR_RETURN ("failed to allocate space for rng state",
			GSL_ENOMEM, 0);
    };

  r->name = T->name;
  r->max = T->max;
  r->size = T->size;
  r->set = T->set;
  r->get = T->get;
  r->get_double = T->get_double;

  gsl_rng_set (r, gsl_rng_default_seed);	/* seed the generator */

  return r;
}

gsl_rng *
gsl_rng_cpy (gsl_rng * dest, const gsl_rng * src)
{
  if (dest->size != src->size)
    {
      dest->state = realloc (dest->state, src->size);

      if (dest->state == 0)
	{
	  GSL_ERROR_RETURN ("failed to reallocate space for rng state",
			    GSL_ENOMEM, 0);
	}
    }

  dest->name = src->name;
  dest->max = src->max;
  dest->size = src->size;
  dest->set = src->set;
  dest->get = src->get;
  dest->get_double = src->get_double;

  memcpy (dest->state, src->state, src->size);

  return dest;
}

gsl_rng *
gsl_rng_clone (const gsl_rng * q)
{
  gsl_rng *r = (gsl_rng *) malloc (sizeof (gsl_rng));

  if (r == 0)
    {
      GSL_ERROR_RETURN ("failed to allocate space for rng struct",
			GSL_ENOMEM, 0);
    };

  r->state = malloc (q->size);

  if (r->state == 0)
    {
      free (r);		/* exception in constructor, avoid memory leak */

      GSL_ERROR_RETURN ("failed to allocate space for rng state",
			GSL_ENOMEM, 0);
    };

  r->name = q->name;
  r->max = q->max;
  r->size = q->size;
  r->set = q->set;
  r->get = q->get;
  r->get_double = q->get_double;

  memcpy (r->state, q->state, q->size);

  return r;
}

void
gsl_rng_set (const gsl_rng * r, unsigned long int seed)
{
  (r->set) (r->state, seed);
}

unsigned long int
gsl_rng_get (const gsl_rng * r)
{
  return (r->get) (r->state);
}

double
gsl_rng_uniform (const gsl_rng * r)
{
  return (r->get_double) (r->state);
}

double
gsl_rng_uniform_pos (const gsl_rng * r)
{
  double x ;
  do
    {
      x = (r->get_double) (r->state) ;
    }
  while (x == 0) ;

  return x ;
}

unsigned long int
gsl_rng_uniform_int (const gsl_rng * r, unsigned long int n)
{
  unsigned long int offset = r->min;
  unsigned long int range = r->max - offset;
  unsigned long int scale = range / n;
  unsigned long int k;

  if (n > range) /*FIXME: test this and move it into gsl_rng.h */
    {
      GSL_ERROR_RETURN ("n exceeds maximum value of generator",
			GSL_EINVAL, 0) ;
    }

  do
    {
      k = (((r->get) (r->state)) - offset) / scale;
    }
  while (k > n);

  return k;
}

unsigned long int
gsl_rng_max (const gsl_rng * r)
{
  return r->max;
}

unsigned long int
gsl_rng_min (const gsl_rng * r)
{
  return r->min;
}

const char *
gsl_rng_name (const gsl_rng * r)
{
  return r->name;
}

void
gsl_rng_print_state (const gsl_rng * r)
{
  size_t i;
  unsigned char *p = (unsigned char *) (r->state);
  const size_t n = r->size;

  for (i = 0; i < n; i++)
    {
      /* FIXME: we're assuming that a char is 8 bits */
      printf ("%.2x", *(p + i));
    }

}

void
gsl_rng_free (gsl_rng * r)
{
  free (r->state);
  free (r);
}
