TYPE (gsl_block) *
FUNCTION (gsl_block, alloc) (const size_t n)
{
  TYPE (gsl_block) * b;

  if (n == 0)
    {
      GSL_ERROR_RETURN ("block length n must be positive integer",
			GSL_EDOM, 0);
    }

  b = (TYPE (gsl_block) *) malloc (sizeof (TYPE (gsl_block)));

  if (b == 0)
    {
      GSL_ERROR_RETURN ("failed to allocate space for block struct",
			GSL_ENOMEM, 0);
    }

  b->data = (ATOMIC *) malloc (MULTIPLICITY * n * sizeof (ATOMIC));

  if (b->data == 0)
    {
      free (b);		/* exception in constructor, avoid memory leak */

      GSL_ERROR_RETURN ("failed to allocate space for block data",
			GSL_ENOMEM, 0);
    }

  b->size = n;

  return b;
}

TYPE (gsl_block) *
FUNCTION (gsl_block, calloc) (const size_t n)
{
  size_t i;

  TYPE (gsl_block) * b = FUNCTION (gsl_block, alloc) (n);

  if (b == 0)
    return 0;

  /* initialize block to zero */

  for (i = 0; i < MULTIPLICITY * n; i++)
    {
      b->data[i] = 0;
    }

  return b;
}

void
FUNCTION (gsl_block, free) (TYPE (gsl_block) * b)
{
  free (b->data);
  free (b);
}
