#include <config.h>

#include <gsl_complex.h>
#include <gsl_errno.h>
#include <gsl_fft_real.h>

#include "fft.h"

int
gsl_fft_real_unpack (const double real_coefficient[],
		     double complex_coefficient[],
		     const size_t n)
{
  size_t i;

  if (n == 0)
    {
      GSL_ERROR ("length n must be positive integer", GSL_EDOM);
    }

  for (i = 0; i < n; i++)
    {
      REAL(complex_coefficient,1,i) = real_coefficient[i];
      IMAG(complex_coefficient,1,i) = 0.0;
    }

  return 0;

}
