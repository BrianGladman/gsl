#include <config.h>
#include <math.h>

#include <gsl_complex.h>
#include <gsl_dft_complex.h>
#include <gsl_errno.h>

int
gsl_dft_complex_forward (const gsl_complex data[],
			 gsl_complex result[],
			 const size_t n)
{
  gsl_fft_direction sign = forward;
  int status = gsl_dft_complex (data, result, n, sign);
  return status;
}

int
gsl_dft_complex_backward (const gsl_complex data[],
			  gsl_complex result[],
			  const size_t n)
{
  gsl_fft_direction sign = backward;
  int status = gsl_dft_complex (data, result, n, sign);
  return status;
}


int
gsl_dft_complex_inverse (const gsl_complex data[],
			 gsl_complex result[],
			 const size_t n)
{
  gsl_fft_direction sign = backward;
  int status = gsl_dft_complex (data, result, n, sign);

  /* normalize inverse fft with 1/n */

  {
    const double norm = 1.0 / n;
    size_t i;
    for (i = 0; i < n; i++)
      {
	result[i].real *= norm;
	result[i].imag *= norm;
      }
  }
  return status;
}

int
gsl_dft_complex (const double data[], double result[],
		 const size_t n, const gsl_fft_direction sign)
{

  size_t i, j, exponent;
  const double d_theta = 2.0 * ((int) sign) * M_PI / (double) n;

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

	  double data_real = ;
	  double data_imag = ;

	  sum_real += w_real * data_real - w_imag * data_imag;
	  sum_imag += w_real * data_imag + w_imag * data_real;

	  exponent = (exponent + i) % n;
	}
      result[i].real = sum_real;
      result[i].imag = sum_imag;
    }
  return 0;
}
