#include "complex_internal.h"

int
FUNCTION(gsl_fft_real,unpack) (const BASE real_coefficient[],
			       BASE complex_coefficient[],
			       const size_t stride, const size_t n)
{
  size_t i;

  if (n == 0)
    {
      GSL_ERROR ("length n must be positive integer", GSL_EDOM);
    }

  for (i = 0; i < n; i++)
    {
      REAL(complex_coefficient,stride,i) = real_coefficient[i * stride];
      IMAG(complex_coefficient,stride,i) = 0.0;
    }

  return 0;

}
