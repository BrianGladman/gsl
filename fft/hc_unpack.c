#include <config.h>
#include <gsl_complex.h>
#include <gsl_errno.h>
#include <gsl_fft_halfcomplex.h>

int
gsl_fft_halfcomplex_unpack (const double halfcomplex_coefficient[],
			    gsl_complex complex_coefficient[],
			    const size_t n)
{
  size_t i;

  if (n == 0)
    {
      GSL_ERROR ("length n must be positive integer", GSL_EDOM);
    }

  GSL_REAL(complex_coefficient[0]) = halfcomplex_coefficient[0];
  GSL_IMAG(complex_coefficient[0]) = 0.0;

  for (i = 1; i < n - i; i++)
    {
      const double hc_real = halfcomplex_coefficient[2 * i - 1];
      const double hc_imag = halfcomplex_coefficient[2 * i];

      GSL_REAL(complex_coefficient[i]) = hc_real;
      GSL_IMAG(complex_coefficient[i]) = hc_imag;
      GSL_REAL(complex_coefficient[n - i]) = hc_real;
      GSL_IMAG(complex_coefficient[n - i]) = -hc_imag;
    }

  if (i == n - i)
    {
      GSL_REAL(complex_coefficient[i]) = halfcomplex_coefficient[n - 1];
      GSL_IMAG(complex_coefficient[i]) = 0.0;
    }

  return 0;

}
