#include <config.h>

#include <gsl_complex.h>
#include <gsl_errno.h>
#include <gsl_fft_real.h>

int
gsl_fft_real_unpack (const double real_coefficient[],
		     gsl_complex complex_coefficient[],
		     const size_t n)
{
  size_t i;

  if (n == 0)
    {
      GSL_ERROR ("length n must be positive integer", GSL_EDOM);
    }

  for (i = 0; i < n; i++)
    {
      complex_coefficient[i].real = real_coefficient[i];
      complex_coefficient[i].imag = 0.0;
    }

  return 0;

}
