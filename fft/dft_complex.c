#include <config.h>
#include <math.h>

#include <gsl_complex.h>
#include <gsl_dft_complex.h>
#include <gsl_errno.h>

#include "fft.h"

int
gsl_dft_complex_forward (const gsl_vector_complex * data,
			 gsl_vector_complex * result)
{
  gsl_fft_direction sign = forward;
  int status = gsl_dft_complex (data, result, sign);
  return status;
}

int
gsl_dft_complex_backward (const gsl_vector_complex * data,
			  gsl_vector_complex * result)
{
  gsl_fft_direction sign = backward;
  int status = gsl_dft_complex (data, result, sign);
  return status;
}


int
gsl_dft_complex_inverse (const gsl_vector_complex * data,
			 gsl_vector_complex * result)
{
  gsl_fft_direction sign = backward;
  int status = gsl_dft_complex (data, result, sign);

  const size_t n = data->size ;
  const size_t ostride = result->stride ;

  double * const out = result->data ;

  /* normalize inverse fft with 1/n */

  {
    const double norm = 1.0 / n;
    size_t i;
    for (i = 0; i < n; i++)
      {
	REAL(out,ostride,i) *= norm;
	IMAG(out,ostride,i) *= norm;
      }
  }
  return status;
}

int
gsl_dft_complex (const gsl_vector_complex * data, 
		 gsl_vector_complex * result,
		 const gsl_fft_direction sign)
{

  size_t i, j, exponent;

  size_t n = data->size ;

  double * const in = data->data ;
  double * const out = result->data ;

  size_t istride = data->stride ;
  size_t ostride = result->stride ;

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

	  double data_real = REAL(in,istride,j);
	  double data_imag = IMAG(in,istride,j);

	  sum_real += w_real * data_real - w_imag * data_imag;
	  sum_imag += w_real * data_imag + w_imag * data_real;

	  exponent = (exponent + i) % n;
	}
      REAL(out,ostride,i) = sum_real;
      IMAG(out,ostride,i) = sum_imag;
    }

  return 0;
}
