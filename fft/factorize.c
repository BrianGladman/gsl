#include <config.h>
#include <gsl_errno.h>
#include <gsl_fft_complex.h>

#include "factorize.h"

int
gsl_fft_complex_factorize (const size_t n,
			   size_t *nf,
			   size_t factors[])
{
  const size_t complex_subtransforms[] =
  {7, 6, 5, 4, 3, 2, 0};

  /* other factors can be added here if their transform modules are
     implemented. The end of the list is marked by 0. */

  int status = gsl_fft_factorize (n, complex_subtransforms, nf, factors);
  return status;
}

int
gsl_fft_halfcomplex_factorize (const size_t n,
			       size_t *nf,
			       size_t factors[])
{
  const size_t halfcomplex_subtransforms[] =
  {5, 4, 3, 2, 0};

  int status = gsl_fft_factorize (n, halfcomplex_subtransforms, nf, factors);
  return status;
}

int
gsl_fft_real_factorize (const size_t n,
			size_t *nf,
			size_t factors[])
{
  const size_t real_subtransforms[] =
  {5, 4, 3, 2, 0};

  int status = gsl_fft_factorize (n, real_subtransforms, nf, factors);
  return status;
}


int
gsl_fft_factorize (const size_t n,
		   const size_t implemented_subtransforms[],
		   size_t *n_factors,
		   size_t factors[])

{
  size_t nf = 0;
  size_t ntest = n;
  size_t factor;
  size_t i = 0;

  if (n == 0)
    {
      GSL_ERROR ("length n must be positive integer", GSL_EDOM);
    }

  if (n == 1)
    {
      factors[0] = 1;
      *n_factors = 1;
      return 0;
    }

  /* deal with the implemented factors first */

  while (implemented_subtransforms[i] && ntest != 1)
    {
      factor = implemented_subtransforms[i];
      while ((ntest % factor) == 0)
	{
	  ntest = ntest / factor;
	  factors[nf] = factor;
	  nf++;
	}
      i++;
    }

  /* deal with any other even prime factors (there is only one) */

  factor = 2;

  while ((ntest % factor) == 0 && (ntest != 1))
    {
      ntest = ntest / factor;
      factors[nf] = factor;
      nf++;
    }

  /* deal with any other odd prime factors */

  factor = 3;

  while (ntest != 1)
    {
      while ((ntest % factor) != 0)
	{
	  factor += 2;
	}
      ntest = ntest / factor;
      factors[nf] = factor;
      nf++;
    }

  /* check that the factorization is correct */
  {
    size_t product = 1;

    for (i = 0; i < nf; i++)
      {
	product *= factors[i];
      }

    if (product != n)
      {
	GSL_ERROR ("factorization failed", GSL_ESANITY);
      }
  }

  *n_factors = nf;

  return 0;
}


int gsl_fft_binary_logn (const size_t n)
{
  size_t ntest ;
  size_t binary_logn = 0 ;
  size_t k = 1;

  while (k < n)
    {
      k *= 2;
      binary_logn++;
    }

  ntest = (1 << binary_logn) ;

  if (n != ntest )       
    {
      return -1 ; /* n is not a power of 2 */
    } 

  return binary_logn;
}




