#include <config.h>
#include <math.h>
#include <stdlib.h>
#include <gsl_math.h>
#include <gsl_complex.h>
#include <gsl_errno.h>

#include <gsl_dft_complex.h>

#include "complex_internal.h"
#include "signals.h"

int
fft_signal_complex_pulse (const size_t k,
			      const size_t n,
			      const double z_real,
			      const double z_imag,
			      double data[],
			      double fft[])
{
  size_t j;

  if (n == 0)
    {
      GSL_ERROR ("length n must be positive integer", GSL_EDOM);
    }

  /* gsl_complex pulse at position k,  data[j] = z * delta_{jk} */

  for (j = 0; j < n; j++)
    {
      REAL(data,1,j) = 0.0;
      IMAG(data,1,j) = 0.0;
    }

  REAL(data,1,k % n) = z_real;
  IMAG(data,1,k % n) = z_imag;

  /* fourier transform, fft[j] = z * exp(-2 pi i j k / n) */

  for (j = 0; j < n; j++)
    {
      const double arg = -2 * M_PI * ((double) ((j * k) % n)) / ((double) n);
      const double w_real = cos (arg);
      const double w_imag = sin (arg);
      REAL(fft,1,j) = w_real * z_real - w_imag * z_imag;
      IMAG(fft,1,j) = w_real * z_imag + w_imag * z_real;
    }

  return 0;

}


int
fft_signal_complex_constant (const size_t n,
				 const double z_real,
				 const double z_imag,
				 double data[],
				 double fft[])
{
  size_t j;

  if (n == 0)
    {
      GSL_ERROR ("length n must be positive integer", GSL_EDOM);
    }

  /* constant, data[j] = z */

  for (j = 0; j < n; j++)
    {
      REAL(data,1,j) = z_real;
      IMAG(data,1,j) = z_imag;
    }

  /* fourier transform, fft[j] = n z delta_{j0} */

  for (j = 0; j < n; j++)
    {
      REAL(fft,1,j) = 0.0;
      IMAG(fft,1,j) = 0.0;
    }

  REAL(fft,1,0) = ((double) n) * z_real;
  IMAG(fft,1,0) = ((double) n) * z_imag;

  return 0;

}


int
fft_signal_complex_exp (const int k,
			    const size_t n,
			    const double z_real,
			    const double z_imag,
			    double data[],
			    double fft[])
{
  size_t j;

  if (n == 0)
    {
      GSL_ERROR ("length n must be positive integer", GSL_EDOM);
    }

  /* exponential,  data[j] = z * exp(2 pi i j k) */

  for (j = 0; j < n; j++)
    {
      const double arg = 2 * M_PI * ((double) ((j * k) % n)) / ((double) n);
      const double w_real = cos (arg);
      const double w_imag = sin (arg);
      REAL(data,1,j) = w_real * z_real - w_imag * z_imag;
      IMAG(data,1,j) = w_real * z_imag + w_imag * z_real;
    }

  /* fourier transform, fft[j] = z * delta{(j - k),0} */

  for (j = 0; j < n; j++)
    {
      REAL(fft,1,j) = 0.0;
      IMAG(fft,1,j) = 0.0;
    }

  {
    int freq;

    if (k <= 0)
      {
	freq = (n-k) % n ;
      }
    else
      {
	freq = (k % n);
      };

    REAL(fft,1,freq) = ((double) n) * z_real;
    IMAG(fft,1,freq) = ((double) n) * z_imag;
  }

  return 0;

}


int
fft_signal_complex_exppair (const int k1,
				     const int k2,
				     const size_t n,
				     const double z1_real,
				     const double z1_imag,
				     const double z2_real,
				     const double z2_imag,
				     double data[],
				     double fft[])
{
  size_t j;

  if (n == 0)
    {
      GSL_ERROR ("length n must be positive integer", GSL_EDOM);
    }

  /* exponential,  data[j] = z1 * exp(2 pi i j k1) + z2 * exp(2 pi i j k2) */

  for (j = 0; j < n; j++)
    {
      const double arg1 = 2 * M_PI * ((double) ((j * k1) % n)) / ((double) n);
      const double w1_real = cos (arg1);
      const double w1_imag = sin (arg1);
      const double arg2 = 2 * M_PI * ((double) ((j * k2) % n)) / ((double) n);
      const double w2_real = cos (arg2);
      const double w2_imag = sin (arg2);
      REAL(data,1,j) = w1_real * z1_real - w1_imag * z1_imag;
      IMAG(data,1,j) = w1_real * z1_imag + w1_imag * z1_real;
      REAL(data,1,j) += w2_real * z2_real - w2_imag * z2_imag;
      IMAG(data,1,j) += w2_real * z2_imag + w2_imag * z2_real;
    }

  /* fourier transform, fft[j] = z1 * delta{(j - k1),0} + z2 *
     delta{(j - k2,0)} */

  for (j = 0; j < n; j++)
    {
      REAL(fft,1,j) = 0.0;
      IMAG(fft,1,j) = 0.0;
    }

  {
    int freq1, freq2;

    if (k1 <= 0)
      {
	freq1 = (n - k1) % n;
      }
    else
      {
	freq1 = (k1 % n);
      };

    if (k2 <= 0)
      {
	freq2 = (n - k2) % n;
      }
    else
      {
	freq2 = (k2 % n);
      };

    REAL(fft,1,freq1) += ((double) n) * z1_real;
    IMAG(fft,1,freq1) += ((double) n) * z1_imag;
    REAL(fft,1,freq2) += ((double) n) * z2_real;
    IMAG(fft,1,freq2) += ((double) n) * z2_imag;
  }

  return 0;

}


int
fft_signal_complex_noise (const size_t n,
			      double data[],
			      double fft[])
{
  size_t i;
  int status;

  if (n == 0)
    {
      GSL_ERROR ("length n must be positive integer", GSL_EDOM);
    }

  for (i = 0; i < n; i++)
    {
      REAL(data,1,i) = ((double) rand ()) / RAND_MAX;
      IMAG(data,1,i) = ((double) rand ()) / RAND_MAX;
    }

  /* compute the dft */
  status = gsl_dft_complex_forward (data, 1, n, fft);

  return status;
}


int
fft_signal_real_noise (const size_t n,
			   double data[],
			   double fft[])
{
  size_t i;
  int status;

  if (n == 0)
    {
      GSL_ERROR ("length n must be positive integer", GSL_EDOM);
    }

  for (i = 0; i < n; i++)
    {
      REAL(data,1,i) = ((double) rand ()) / RAND_MAX;
      IMAG(data,1,i) = 0.0;
    }

  /* compute the dft */
  status = gsl_dft_complex_forward (data, 1, n, fft);

  return status;
}

