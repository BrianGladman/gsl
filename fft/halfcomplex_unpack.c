#include <gsl_complex.h>
#include <gsl_errno.h>
#include <gsl_fft_halfcomplex.h>

int
gsl_fft_halfcomplex_unpack (const double halfcomplex_coefficient[],
			    gsl_complex complex_coefficient[],
			    const unsigned int n)
{
  unsigned int i;

  if (n == 0)
    {
      GSL_ERROR ("length n must be positive integer", GSL_EDOM);
    }

  complex_coefficient[0].real = halfcomplex_coefficient[0];
  complex_coefficient[0].imag = 0.0;

  for (i = 1; i < n - i; i++)
    {
      const double hc_real = halfcomplex_coefficient[2 * i - 1];
      const double hc_imag = halfcomplex_coefficient[2 * i];

      complex_coefficient[i].real = hc_real;
      complex_coefficient[i].imag = hc_imag;
      complex_coefficient[n - i].real = hc_real;
      complex_coefficient[n - i].imag = -hc_imag;
    }

  if (i == n - i)
    {
      complex_coefficient[i].real = halfcomplex_coefficient[n - 1];
      complex_coefficient[i].imag = 0.0;
    }

  return 0;

}
