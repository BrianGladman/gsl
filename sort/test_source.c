void TYPE (test_sort_vector) (size_t N, size_t stride);
void FUNCTION (my, initialize) (TYPE (gsl_vector) * v);
void FUNCTION (my, randomize) (TYPE (gsl_vector) * v);
int FUNCTION (my, check) (TYPE (gsl_vector) * data, TYPE (gsl_vector) * orig);
int FUNCTION (my, pcheck) (gsl_permutation * p, TYPE (gsl_vector) * data, TYPE (gsl_vector) * orig);

void
TYPE (test_sort_vector) (size_t N, size_t stride)
{
  int status;

  TYPE (gsl_block) * b1 = FUNCTION (gsl_block, calloc) (N * stride);
  TYPE (gsl_block) * b2 = FUNCTION (gsl_block, calloc) (N * stride);

  TYPE (gsl_vector) * orig = FUNCTION (gsl_vector, alloc_from_block) (b1, 0, N, stride);
  TYPE (gsl_vector) * data = FUNCTION (gsl_vector, alloc_from_block) (b2, 0, N, stride);

  gsl_permutation *p = gsl_permutation_alloc (N);

  FUNCTION (my, initialize) (orig);

  /* Already sorted */
  FUNCTION (gsl_vector, memcpy) (data, orig);

  status = FUNCTION (gsl_sort_vector, index) (p, data);
  status |= FUNCTION (my, pcheck) (p, data, orig);
  gsl_test (status, "indexing " NAME (gsl_vector) ", n = %u, stride = %u, ordered", N, stride);

  TYPE (gsl_sort_vector) (data);
  status = FUNCTION (my, check) (data, orig);
  gsl_test (status, "sorting, " NAME (gsl_vector) ", n = %u, stride = %u, ordered", N, stride);

  /* Reverse the data */

  FUNCTION (gsl_vector, memcpy) (data, orig);
  FUNCTION (gsl_vector, reverse) (data);

  status = FUNCTION (gsl_sort_vector, index) (p, data);
  status |= FUNCTION (my, pcheck) (p, data, orig);
  gsl_test (status, "indexing " NAME (gsl_vector) ", n = %u, stride = %u, reversed", N, stride);

  TYPE (gsl_sort_vector) (data);
  status = FUNCTION (my, check) (data, orig);
  gsl_test (status, "sorting, " NAME (gsl_vector) ", n = %u, stride = %u, reversed", N, stride);

  /* Perform some shuffling */

  FUNCTION (gsl_vector, memcpy) (data, orig);
  FUNCTION (my, randomize) (data);

  status = FUNCTION (gsl_sort_vector, index) (p, data);
  status |= FUNCTION (my, pcheck) (p, data, orig);
  gsl_test (status, "indexing " NAME (gsl_vector) ", n = %u, stride = %u, randomized", N, stride);

  TYPE (gsl_sort_vector) (data);
  status = FUNCTION (my, check) (data, orig);
  gsl_test (status, "sorting, " NAME (gsl_vector) ", n = %u, stride = %u, randomized", N, stride);

  FUNCTION (gsl_vector, free) (orig);
  FUNCTION (gsl_vector, free) (data);
  FUNCTION (gsl_block, free) (b1);
  FUNCTION (gsl_block, free) (b2);
  gsl_permutation_free (p);
}


void
FUNCTION (my, initialize) (TYPE (gsl_vector) * v)
{
  size_t i;
  ATOMIC k = 0, kk;

  /* Must be sorted initially */

  for (i = 0; i < v->size; i++)
    {
      kk = k;
      k++;
      if (k < kk)		/* prevent overflow */
	k = kk;
      FUNCTION (gsl_vector, set) (v, i, k);
    }
}

void
FUNCTION (my, randomize) (TYPE (gsl_vector) * v)
{
  size_t i;

  for (i = 0; i < v->size; i++)
    {
      size_t j = urand () * v->size;
      FUNCTION (gsl_vector, swap_elements) (v, i, j);
    }
}

int
FUNCTION (my, check) (TYPE (gsl_vector) * data, TYPE (gsl_vector) * orig)
{
  size_t i;

  for (i = 0; i < data->size; i++)
    {
      if (FUNCTION (gsl_vector, get) (data, i) != FUNCTION (gsl_vector, get) (orig, i))
	{
	  return GSL_FAILURE;
	}
    }

  return GSL_SUCCESS;
}

int
FUNCTION (my, pcheck) (gsl_permutation * p, TYPE (gsl_vector) * data, TYPE (gsl_vector) * orig)
{
  size_t i;

  for (i = 0; i < p->size; i++)
    {
      if (FUNCTION (gsl_vector, get) (data, p->data[i]) != FUNCTION (gsl_vector, get) (orig, i))
	{
	  return GSL_FAILURE;
	}
    }

  return GSL_SUCCESS;
}


