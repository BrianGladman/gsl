/* Randomly permute (shuffle) N indices

   Supply an integer array x[N], and on return, it will be filled with
   indices 0...N-1 in random order.  The algorithm is from Knuth,
   SemiNumerical Algorithms, v2, p139 */

#include <stdlib.h>
#include <stdio.h>		/* defines NULL */
#include <math.h>		/* defines floor() */
#include <gsl_rng.h>
#include "gsl_randist.h"

int *
gsl_ran_shuffle (const gsl_rng * r, int N, int *x)
{
  int i, k, tmp;

  /* First, do a bunch of memory allocation stuff */
  if (N < 0)
    return NULL;
  if (x == NULL && N > 0)
    {
      x = (int *) calloc ((size_t) N, sizeof (int));
      if (x == NULL)
	return NULL;
      for (i = 0; i < N; ++i)
	x[i] = i;
    }
  if (x != NULL && N == 0)
    {
      free ((char *) x);
      return NULL;
    }
  /* Now here's the algorithm, more or less transcribed
   * from Knuth, who cites Moses and Oakford, and Durstenfeld */

  for (i = N - 1; i >= 0; --i)
    {
      k = floor (i * gsl_rng_uniform (r));
      tmp = x[k];
      x[k] = x[i];
      x[i] = tmp;
    }
  return x;
}
int *
gsl_ran_choose (const gsl_rng * r, int K, int N, int *x)
{
  int n, k;
  /* Choose K out of N items */
  /* return an array x[] of the indices of the N items */
  /* these items will be in sorted order -- you can use
   * shuffle() to randomize them if you wish */

  /* First, do a bunch of memory allocation stuff */
  if (N < K || K < 0)
    {
      if (x != NULL)
	free ((char *) x);
      return NULL;
    }
  if (x == NULL && K > 0)
    {
      x = (int *) calloc ((size_t) K, sizeof (int));
      if (x == NULL)
	return NULL;
    }
  /* Here is the guts of the algorithm: three lines!! */
  for (n = 0, k = 0; n < N && k < K; ++n)
    {
      if ((N - n) * gsl_rng_uniform (r) < K - k)
	{
	  x[k++] = n;
	}
    }
  return x;
}
