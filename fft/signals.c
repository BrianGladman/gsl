#include <config.h>
#include <math.h>
#include <stdlib.h>
#include <gsl_math.h>
#include <gsl_complex.h>
#include <gsl_errno.h>

#include <gsl_fft_signals.h>
#include <gsl_dft_complex.h>

int
gsl_fft_signal_complex_pulse (const unsigned int k,
			      const unsigned int n,
			      const double z_real,
			      const double z_imag,
			      gsl_complex data[],
			      gsl_complex fft[])
{
  unsigned int j;

  if (n == 0)
    {
      GSL_ERROR ("length n must be positive integer", GSL_EDOM);
    }

  /* gsl_complex pulse at position k,  data[j] = z * delta_{jk} */

  for (j = 0; j < n; j++)
    {
      data[j].real = 0.0;
      data[j].imag = 0.0;
    }

  data[k % n].real = z_real;
  data[k % n].imag = z_imag;

  /* fourier transform, fft[j] = z * exp(-2 pi i j k / n) */

  for (j = 0; j < n; j++)
    {
      const double arg = -2 * M_PI * ((double) ((j * k) % n)) / ((double) n);
      const double w_real = cos (arg);
      const double w_imag = sin (arg);
      fft[j].real = w_real * z_real - w_imag * z_imag;
      fft[j].imag = w_real * z_imag + w_imag * z_real;
    }

  return 0;

}


int
gsl_fft_signal_complex_constant (const unsigned int n,
				      const double z_real,
				      const double z_imag,
				      gsl_complex data[],
				      gsl_complex fft[])
{
  unsigned int j;

  if (n == 0)
    {
      GSL_ERROR ("length n must be positive integer", GSL_EDOM);
    }

  /* constant, data[j] = z */

  for (j = 0; j < n; j++)
    {
      data[j].real = z_real;
      data[j].imag = z_imag;
    }

  /* fourier transform, fft[j] = n z delta_{j0} */

  for (j = 0; j < n; j++)
    {
      fft[j].real = 0.0;
      fft[j].imag = 0.0;
    }

  fft[0].real = ((double) n) * z_real;
  fft[0].imag = ((double) n) * z_imag;

  return 0;

}


int
gsl_fft_signal_complex_exp (const int k,
				 const unsigned int n,
				 const double z_real,
				 const double z_imag,
				 gsl_complex data[],
				 gsl_complex fft[])
{
  unsigned int j;

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
      data[j].real = w_real * z_real - w_imag * z_imag;
      data[j].imag = w_real * z_imag + w_imag * z_real;
    }

  /* fourier transform, fft[j] = z * delta{(j - k),0} */

  for (j = 0; j < n; j++)
    {
      fft[j].real = 0.0;
      fft[j].imag = 0.0;
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

    fft[freq].real = ((double) n) * z_real;
    fft[freq].imag = ((double) n) * z_imag;
  }

  return 0;

}


int
gsl_fft_signal_complex_exppair (const int k1,
				     const int k2,
				     const unsigned int n,
				     const double z1_real,
				     const double z1_imag,
				     const double z2_real,
				     const double z2_imag,
				     gsl_complex data[],
				     gsl_complex fft[])
{
  unsigned int j;

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
      data[j].real = w1_real * z1_real - w1_imag * z1_imag;
      data[j].imag = w1_real * z1_imag + w1_imag * z1_real;
      data[j].real += w2_real * z2_real - w2_imag * z2_imag;
      data[j].imag += w2_real * z2_imag + w2_imag * z2_real;
    }

  /* fourier transform, fft[j] = z1 * delta{(j - k1),0} + z2 *
     delta{(j - k2,0)} */

  for (j = 0; j < n; j++)
    {
      fft[j].real = 0.0;
      fft[j].imag = 0.0;
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

    fft[freq1].real += ((double) n) * z1_real;
    fft[freq1].imag += ((double) n) * z1_imag;
    fft[freq2].real += ((double) n) * z2_real;
    fft[freq2].imag += ((double) n) * z2_imag;
  }

  return 0;

}


int
gsl_fft_signal_complex_noise (const unsigned int n,
				   gsl_complex data[],
				   gsl_complex fft[])
{
  unsigned int i;
  int status;

  if (n == 0)
    {
      GSL_ERROR ("length n must be positive integer", GSL_EDOM);
    }

  for (i = 0; i < n; i++)
    {
      data[i].real = ((double) rand ()) / RAND_MAX;
      data[i].imag = ((double) rand ()) / RAND_MAX;
    }

  /* compute the dft */
  status = gsl_dft_complex_forward (data, fft, n);

  return status;
}


int
gsl_fft_signal_real_noise (const unsigned int n,
				gsl_complex data[],
				gsl_complex fft[])
{
  unsigned int i;
  int status;

  if (n == 0)
    {
      GSL_ERROR ("length n must be positive integer", GSL_EDOM);
    }

  for (i = 0; i < n; i++)
    {
      data[i].real = ((double) rand ()) / RAND_MAX;
      data[i].imag = 0.0;
    }

  /* compute the dft */
  status = gsl_dft_complex_forward (data, fft, n);

  return status;
}

