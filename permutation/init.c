#include <config.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_permutation.h>

gsl_permutation *
gsl_permutation_alloc (const size_t n)
{
  gsl_permutation * p;

  if (n == 0)
    {
      GSL_ERROR_RETURN ("permutation length n must be positive integer",
			GSL_EDOM, 0);
    }

  p = (gsl_permutation *) malloc (sizeof (gsl_permutation));

  if (p == 0)
    {
      GSL_ERROR_RETURN ("failed to allocate space for permutation struct",
			GSL_ENOMEM, 0);
    }

  p->data = (size_t *) malloc (n * sizeof (size_t));

  if (p->data == 0)
    {
      free (p);		/* exception in constructor, avoid memory leak */

      GSL_ERROR_RETURN ("failed to allocate space for permutation data",
			GSL_ENOMEM, 0);
    }

  p->size = n;

  return p;
}

gsl_permutation *
gsl_permutation_calloc (const size_t n)
{
  size_t i;

  gsl_permutation * p =  gsl_permutation_alloc (n);

  if (p == 0)
    return 0;

  /* initialize permutation to identity */

  for (i = 0; i < n; i++)
    {
      p->data[i] = i;
    }

  return p;
}

void
gsl_permutation_init (gsl_permutation * p)
{
  const size_t n = p->size ;
  size_t i;

  /* initialize permutation to identity */

  for (i = 0; i < n; i++)
    {
      p->data[i] = i;
    }
}

void
gsl_permutation_free (gsl_permutation * p)
{
  free (p->data);
  free (p);
}
