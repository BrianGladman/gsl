TYPE (gsl_block) *
FUNCTION (gsl_block, alloc) (const size_t n)
{
  TYPE (gsl_block) * v;

  if (n == 0)
    {
      GSL_ERROR_RETURN ("block length n must be positive integer",
			GSL_EDOM, 0);
    }

  v = (TYPE (gsl_block) *) malloc (sizeof (TYPE (gsl_block)));

  if (v == 0)
    {
      GSL_ERROR_RETURN ("failed to allocate space for block struct",
			GSL_ENOMEM, 0);
    }

  v->data = (ATOMIC *) malloc (MULTIPLICITY * n * sizeof (ATOMIC));

  if (v->data == 0)
    {
      free (v);		/* exception in constructor, avoid memory leak */

      GSL_ERROR_RETURN ("failed to allocate space for block data",
			GSL_ENOMEM, 0);
    }

  v->size = n;

  return v;
}

TYPE (gsl_block) *
FUNCTION (gsl_block, calloc) (const size_t n)
{
  size_t i;

  TYPE (gsl_block) * v = FUNCTION (gsl_block, alloc) (n);

  if (v == 0)
    return 0;

  /* initialize block to zero */

  for (i = 0; i < MULTIPLICITY * n; i++)
    {
      v->data[i] = 0;
    }

  return v;
}

void
FUNCTION (gsl_block, free) (TYPE (gsl_block) * v)
{
  free (v->data);
  free (v);
}
