#include <config.h>
#include <gsl_complex.h>
#include <gsl_errno.h>
#include <gsl_fft_halfcomplex.h>

#include "fft.h"

int
gsl_fft_halfcomplex_unpack (const double halfcomplex_coefficient[],
			    double complex_coefficient[],
			    const size_t n)
{
  size_t i;

  if (n == 0)
    {
      GSL_ERROR ("length n must be positive integer", GSL_EDOM);
    }

  REAL(complex_coefficient,1,0) = halfcomplex_coefficient[0];
  IMAG(complex_coefficient,1,0) = 0.0;

  for (i = 1; i < n - i; i++)
    {
      const double hc_real = halfcomplex_coefficient[2 * i - 1];
      const double hc_imag = halfcomplex_coefficient[2 * i];

      REAL(complex_coefficient,1,i) = hc_real;
      IMAG(complex_coefficient,1,i) = hc_imag;
      REAL(complex_coefficient,1,n - i) = hc_real;
      IMAG(complex_coefficient,1,n - i) = -hc_imag;
    }

  if (i == n - i)
    {
      REAL(complex_coefficient,1,i) = halfcomplex_coefficient[n - 1];
      IMAG(complex_coefficient,1,i) = 0.0;
    }

  return 0;

}

