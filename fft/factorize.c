#include <gsl_errno.h>
#include <gsl_fft_complex.h>

int
gsl_fft_complex_factorize (const unsigned int n,
			   unsigned int *nf,
			   unsigned int factors[])
{
  const unsigned int complex_subtransforms[] =
  {6, 5, 4, 3, 2, 0};

  /* other factors can be added here if their transform modules are
     implemented. The end of the list is marked by 0. */

  int status = gsl_fft_factorize (n, complex_subtransforms, nf, factors);
  return status;
}

int
gsl_fft_halfcomplex_factorize (const unsigned int n,
			       unsigned int *nf,
			       unsigned int factors[])
{
  const unsigned int halfcomplex_subtransforms[] =
  {5, 4, 3, 2, 0};

  int status = gsl_fft_factorize (n, halfcomplex_subtransforms, nf, factors);
  return status;
}

int
gsl_fft_real_factorize (const unsigned int n,
			unsigned int *nf,
			unsigned int factors[])
{
  const unsigned int real_subtransforms[] =
  {5, 4, 3, 2, 0};

  int status = gsl_fft_factorize (n, real_subtransforms, nf, factors);
  return status;
}


int
gsl_fft_factorize (const unsigned int n,
		   const unsigned int implemented_subtransforms[],
		   unsigned int *n_factors,
		   unsigned int factors[])

{
  unsigned int nf = 0;
  unsigned int ntest = n;
  unsigned int factor_sum = 0;
  unsigned int factor;
  unsigned int i = 0;

  if (n == 0)
    {
      GSL_ERROR ("length n must be positive integer", GSL_EDOM);
    }

  if (n == 1)
    {
      factors[0] = 1;
      *n_factors = 1;
      return 1;
    }

  /* deal with the implemented factors first */

  while (implemented_subtransforms[i] && ntest != 1)
    {
      factor = implemented_subtransforms[i];
      while ((ntest % factor) == 0)
	{
	  ntest = ntest / factor;
	  factors[nf] = factor;
	  factor_sum += factor;
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
      factor_sum += factor;
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
      factor_sum += factor;
      nf++;
    }

  /* check that the factorization is correct */
  {
    unsigned int product = 1;

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

  /* the sum of the factors gives the scaling of the algorithm

     T ~ O(n factor_sum)

     a well factorized length has a factor sum which is much less than n */

  return factor_sum;
}


int gsl_fft_binary_logn (const unsigned int n)
{
  unsigned int ntest ;
  unsigned int binary_logn = 0 ;
  unsigned int k = 1;

  while (k < n)
    {
      k *= 2;
      binary_logn++;
    }

  ntest = (1 << binary_logn) ;

  if (n != ntest )
    {
      /* n is not a power of 2 */
      return -1 ; 
    } 
  else 
    {
      return binary_logn;
    }
      
}




