#include <math.h>
#include <gsl_complex.h>
#include <gsl_dft_complex.h>
#include <gsl_errno.h>

int
gsl_dft_complex_forward (const complex data[],
			 complex result[],
			 const unsigned int n)
{
  gsl_fft_direction sign = forward;
  int status = gsl_dft_complex (data, result, n, sign);
  return status;
}

int
gsl_dft_complex_backward (const complex data[],
			  complex result[],
			  const unsigned int n)
{
  gsl_fft_direction sign = backward;
  int status = gsl_dft_complex (data, result, n, sign);
  return status;
}


int
gsl_dft_complex_inverse (const complex data[],
			 complex result[],
			 const unsigned int n)
{
  gsl_fft_direction sign = backward;
  int status = gsl_dft_complex (data, result, n, sign);

  /* normalize inverse fft with 1/n */

  {
    const double norm = 1.0 / n;
    unsigned int i;
    for (i = 0; i < n; i++)
      {
	result[i].real *= norm;
	result[i].imag *= norm;
      }
  }
  return status;
}

int
gsl_dft_complex (const complex data[],
		 complex result[],
		 const unsigned int n,
		 const gsl_fft_direction sign)
{

  unsigned int i, j, exponent;
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

	  double data_real = data[j].real;
	  double data_imag = data[j].imag;

	  sum_real += w_real * data_real - w_imag * data_imag;
	  sum_imag += w_real * data_imag + w_imag * data_real;

	  exponent = (exponent + i) % n;
	}
      result[i].real = sum_real;
      result[i].imag = sum_imag;
    }
  return 0;
}
