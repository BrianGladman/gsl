#include <config.h>
#include <gsl/gsl_fft.h>

#include "complex_internal.h"
#include "bitreverse.h"

static int 
FUNCTION(fft_complex,bitreverse_order) (BASE data[], 
					const size_t stride,
					const size_t n,
					size_t logn)
{
  /* This is the Goldrader bit-reversal algorithm */

  size_t i;
  size_t j = 0;

  logn = 0 ; /* not needed for this algorithm */

  for (i = 0; i < n - 1; i++)
    {
      size_t k = n / 2 ;

      if (i < j)
	{
	  const BASE tmp_real = REAL(data,stride,i);
	  const BASE tmp_imag = IMAG(data,stride,i);
	  REAL(data,stride,i) = REAL(data,stride,j);
	  IMAG(data,stride,i) = IMAG(data,stride,j);
	  REAL(data,stride,j) = tmp_real;
	  IMAG(data,stride,j) = tmp_imag;
	}

      while (k <= j) 
	{
	  j = j - k ;
	  k = k / 2 ;
	}

      j += k ;
    }

  return 0;
}


static int 
FUNCTION(fft_real,bitreverse_order) (BASE data[], 
				const size_t stride, 
				const size_t n,
				size_t logn)
{
  /* This is the Goldrader bit-reversal algorithm */

  size_t i;
  size_t j = 0;

  logn = 0 ; /* not needed for this algorithm */

  for (i = 0; i < n - 1; i++)
    {
      size_t k = n / 2 ;

      if (i < j)
	{
	  const BASE tmp = VECTOR(data,stride,i);
	  VECTOR(data,stride,i) = VECTOR(data,stride,j);
	  VECTOR(data,stride,j) = tmp;
	}

      while (k <= j) 
	{
	  j = j - k ;
	  k = k / 2 ;
	}

      j += k ;
    }

  return 0;
}

