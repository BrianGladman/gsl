#include <stdlib.h>
#include <math.h>

#include <gsl_errno.h>
#include <gsl_complex.h>
#include <gsl_fft_halfcomplex.h>

#include <fft_halfcomplex.h>

int
gsl_fft_halfcomplex_backward (double *data,
			      const unsigned int n,
			      const gsl_fft_halfcomplex_wavetable * wavetable)
{
  int status = gsl_fft_halfcomplex (data, n, wavetable) ;
  return status ;
}

int
gsl_fft_halfcomplex_inverse (double * data,
			     const unsigned int n,
			     const gsl_fft_halfcomplex_wavetable * wavetable)
{
  int status = gsl_fft_halfcomplex (data, n, wavetable);

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
	data[i] *= norm;
      }
  }
  return status;
}

int
gsl_fft_halfcomplex (double *data,
		     const unsigned int n,
		     const gsl_fft_halfcomplex_wavetable * wavetable)
{

  unsigned int factor, product, q, state;
  unsigned int i;
  unsigned int nf;
  int product_1;
  int tskip;
  double *from, *to, *scratch;
  complex *twiddle1, *twiddle2, *twiddle3, *twiddle4;

  if (n == 0)
    {
      GSL_ERROR ("length n must be positive integer", GSL_EDOM);
    }

  if (n == 1)
    {				/* FFT of one data point is the identity */
      return 0;
    }

  if (n != wavetable->n)
    {
      GSL_ERROR ("wavetable does not match length of data", GSL_EINVAL);
    }

  nf = wavetable->nf;

  scratch = (double *) (wavetable->scratch);
  product = 1;

  state = 0;
  from = data;
  to = scratch;

  for (i = 0; i < nf; i++)
    {
      factor = wavetable->factor[i];
      product_1 = product;
      product *= factor;
      q = n / product;

      tskip = (q + 1) / 2 - 1;

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
	  gsl_fft_halfcomplex_pass_2 (from, to, product, n, twiddle1);
	}
      else if (factor == 3)
	{
	  twiddle1 = wavetable->twiddle[i];
	  twiddle2 = twiddle1 + tskip;
	  gsl_fft_halfcomplex_pass_3 (from, to, product, n, twiddle1,
				      twiddle2);
	}
      else if (factor == 4)
	{
	  twiddle1 = wavetable->twiddle[i];
	  twiddle2 = twiddle1 + tskip;
	  twiddle3 = twiddle2 + tskip;
	  gsl_fft_halfcomplex_pass_4 (from, to, product, n, twiddle1,
				      twiddle2, twiddle3);
	}
      else if (factor == 5)
	{
	  twiddle1 = wavetable->twiddle[i];
	  twiddle2 = twiddle1 + tskip;
	  twiddle3 = twiddle2 + tskip;
	  twiddle4 = twiddle3 + tskip;
	  gsl_fft_halfcomplex_pass_5 (from, to, product, n, twiddle1,
				      twiddle2, twiddle3, twiddle4);
	}
      else
	{
	  twiddle1 = wavetable->twiddle[i];
	  gsl_fft_halfcomplex_pass_n (from, to, factor, product, n,
				      twiddle1);
	}
    }

  if (state == 1)		/* copy results back from scratch to data */
    {
      memcpy (data, scratch, n * sizeof (double));
    }

  return 0;

}
