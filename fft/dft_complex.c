#include <config.h>
#include <math.h>

#include <gsl_complex.h>
#include <gsl_dft_complex.h>
#include <gsl_errno.h>

#include "complex_internal.h"

int
gsl_dft_complex_forward (const double data[], 
			 const size_t stride, const size_t n,
			 double result[])
{
  gsl_fft_direction sign = forward;
  int status = gsl_dft_complex (data, stride, n, result, sign);
  return status;
}

int
gsl_dft_complex_backward (const double data[], 
			 const size_t stride, const size_t n,
			 double result[])
{
  gsl_fft_direction sign = backward;
  int status = gsl_dft_complex (data, stride, n, result, sign);
  return status;
}


int
gsl_dft_complex_inverse (const double data[], 
			 const size_t stride, const size_t n,
			 double result[])
{
  gsl_fft_direction sign = backward;
  int status = gsl_dft_complex (data, stride, n, result, sign);

  /* normalize inverse fft with 1/n */

  {
    const double norm = 1.0 / n;
    size_t i;
    for (i = 0; i < n; i++)
      {
	REAL(result,1,i) *= norm;
	IMAG(result,1,i) *= norm;
      }
  }
  return status;
}

int
gsl_dft_complex (const double data[], 
		 const size_t stride, const size_t n,
		 double result[],
		 const gsl_fft_direction sign)
{

  size_t i, j, exponent;

  const double d_theta = 2.0 * ((int) sign) * M_PI / (double) n;

  /* FIXME: check that input length == output length and give error */

  for (i = 0; i < n; i++)
    {
      double sum_real = 0;
      double sum_imag = 0;

      exponent = 0;

      for (j = 0; j < n; j++)
	{
	  double theta = d_theta * (double) exponent;
	  /* sum = exp(i theta) * data[j] */

	  double w_real = cos (theta);
	  double w_imag = sin (theta);

	  double data_real = REAL(data,stride,j);
	  double data_imag = IMAG(data,stride,j);

	  sum_real += w_real * data_real - w_imag * data_imag;
	  sum_imag += w_real * data_imag + w_imag * data_real;

	  exponent = (exponent + i) % n;
	}
      REAL(result,1,i) = sum_real;
      IMAG(result,1,i) = sum_imag;
    }

  return 0;
}
