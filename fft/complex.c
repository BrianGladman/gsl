#include <stdlib.h>
#include <math.h>

#include <gsl_errno.h>
#include <gsl_complex.h>
#include <gsl_fft_complex.h>

int
gsl_fft_complex_forward (complex data[],
			 const unsigned int n,
			 const gsl_fft_complex_wavetable * wavetable)
{
  gsl_fft_direction sign = forward;
  int status = gsl_fft_complex (data, n, wavetable, sign);
  return status;
}

int
gsl_fft_complex_backward (complex data[],
			  const unsigned int n,
			  const gsl_fft_complex_wavetable * wavetable)
{
  gsl_fft_direction sign = backward;
  int status = gsl_fft_complex (data, n, wavetable, sign);
  return status;
}

int
gsl_fft_complex_inverse (complex data[],
			 const unsigned int n,
			 const gsl_fft_complex_wavetable * wavetable)
{
  gsl_fft_direction sign = backward;
  int status = gsl_fft_complex (data, n, wavetable, sign);

  if (status)
    {
      return status;
    }

  /* normalize inverse fft with 1/n */

  {
    const double norm = 1.0 / n;
    unsigned int i;
    for (i = 0; i < n; i++)
      {
	data[i].real *= norm;
	data[i].imag *= norm;
      }
  }
  return status;
}

int
gsl_fft_complex (complex data[],
		 const unsigned int n,
		 const gsl_fft_complex_wavetable * wavetable,
		 const gsl_fft_direction sign)
{

  const unsigned int nf = wavetable->nf;

  unsigned int i;

  unsigned int q, product = 1;

  complex *scratch = wavetable->scratch;
  complex *twiddle1, *twiddle2, *twiddle3, *twiddle4, *twiddle5, *twiddle6;

  unsigned int state = 0;
  complex *from = data;
  complex *to = scratch;

  if (n == 0)
    {
      GSL_ERROR ("length n must be positive integer", GSL_EDOM);
    }

  if (n == 1)
    {				/* FFT of 1 data point is the identity */
      return 0;
    }

  if (n != wavetable->n)
    {
      GSL_ERROR ("wavetable does not match length of data", GSL_EINVAL);
    }


  for (i = 0; i < nf; i++)
    {
      const unsigned int factor = wavetable->factor[i];
      product *= factor;
      q = n / product;

      if (state == 0)
	{
	  from = data;
	  to = scratch;
	  state = 1;
	}
      else
	{
	  from = scratch;
	  to = data;
	  state = 0;
	}

      if (factor == 2)
	{
	  twiddle1 = wavetable->twiddle[i];
	  gsl_fft_complex_pass_2 (from, to, sign, product, n, twiddle1);
	}
      else if (factor == 3)
	{
	  twiddle1 = wavetable->twiddle[i];
	  twiddle2 = twiddle1 + q;
	  gsl_fft_complex_pass_3 (from, to, sign, product, n, twiddle1,
				  twiddle2);
	}
      else if (factor == 4)
	{
	  twiddle1 = wavetable->twiddle[i];
	  twiddle2 = twiddle1 + q;
	  twiddle3 = twiddle2 + q;
	  gsl_fft_complex_pass_4 (from, to, sign, product, n, twiddle1,
				  twiddle2, twiddle3);
	}
      else if (factor == 5)
	{
	  twiddle1 = wavetable->twiddle[i];
	  twiddle2 = twiddle1 + q;
	  twiddle3 = twiddle2 + q;
	  twiddle4 = twiddle3 + q;
	  gsl_fft_complex_pass_5 (from, to, sign, product, n, twiddle1,
				  twiddle2, twiddle3, twiddle4);
	}
      else if (factor == 6)
	{
	  twiddle1 = wavetable->twiddle[i];
	  twiddle2 = twiddle1 + q;
	  twiddle3 = twiddle2 + q;
	  twiddle4 = twiddle3 + q;
	  twiddle5 = twiddle4 + q;
	  gsl_fft_complex_pass_6 (from, to, sign, product, n, twiddle1,
				  twiddle2, twiddle3, twiddle4, twiddle5);
	}
      else if (factor == 7)
	{
	  twiddle1 = wavetable->twiddle[i];
	  twiddle2 = twiddle1 + q;
	  twiddle3 = twiddle2 + q;
	  twiddle4 = twiddle3 + q;
	  twiddle5 = twiddle4 + q;
	  twiddle6 = twiddle5 + q;
	  gsl_fft_complex_pass_7 (from, to, sign, product, n, twiddle1,
				  twiddle2, twiddle3, twiddle4, twiddle5, twiddle6);
	}
      else
	{
	  twiddle1 = wavetable->twiddle[i];
	  gsl_fft_complex_pass_n (from, to, sign, factor, product, n,
				  twiddle1);
	}
    }

  if (state == 1)		/* copy results back from scratch to data */
    {
      memcpy (data, scratch, n * sizeof (complex));
    }

  return 0;

}
