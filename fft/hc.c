#include <config.h>
#include <stdlib.h>
#include <math.h>

#include <gsl_errno.h>
#include <gsl_complex.h>
#include <gsl_vector.h>
#include <gsl_fft_halfcomplex.h>

#include "fft_halfcomplex.h"

int
gsl_fft_halfcomplex_backward (gsl_vector * data,
			      const gsl_fft_halfcomplex_wavetable * wavetable)
{
  int status = gsl_fft_halfcomplex (data, wavetable) ;
  return status ;
}

int
gsl_fft_halfcomplex_inverse (gsl_vector * data,
			     const gsl_fft_halfcomplex_wavetable * wavetable)
{
  int status = gsl_fft_halfcomplex (data, wavetable);

  if (status)
    {
      return status;
    }

  /* normalize inverse fft with 1/n */

  {
    const size_t n = data->size ;
    const size_t stride = data->stride ;
    double * const out = data->data ;

    const double norm = 1.0 / n;
    size_t i;
    for (i = 0; i < n; i++)
      {
	out[stride*i] *= norm;
      }
  }
  return status;
}

int
gsl_fft_halfcomplex (gsl_vector * data,
		     const gsl_fft_halfcomplex_wavetable * wavetable)
{
  size_t factor, product, q, state;
  size_t i;
  size_t nf;
  int product_1;
  int tskip;
  gsl_complex *twiddle1, *twiddle2, *twiddle3, *twiddle4;

  const size_t n = data->size ;

  double * const a = data->data;
  double * const b = wavetable->scratch;

  const size_t astride = data->stride ;
  const size_t bstride = 1 ;

  double * in = a;
  size_t istride = astride;
  double * out = b;
  size_t ostride = bstride;

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
  product = 1;
  state = 0;

  for (i = 0; i < nf; i++)
    {
      factor = wavetable->factor[i];
      product_1 = product;
      product *= factor;
      q = n / product;

      tskip = (q + 1) / 2 - 1;

      if (state == 0)
	{
	  in = a;
	  istride = astride;
	  out = b;
	  ostride = bstride;
	  state = 1;
	}
      else
	{
	  in = b;
	  istride = bstride;
	  out = a;
	  ostride = astride;
	  state = 0;
	}

      if (factor == 2)
	{
	  twiddle1 = wavetable->twiddle[i];
	  gsl_fft_halfcomplex_pass_2 (in, istride, out, ostride, 
				      product, n, twiddle1);
	}
      else if (factor == 3)
	{
	  twiddle1 = wavetable->twiddle[i];
	  twiddle2 = twiddle1 + tskip;
	  gsl_fft_halfcomplex_pass_3 (in, istride, out, ostride,
				      product, n, twiddle1, twiddle2);
	}
      else if (factor == 4)
	{
	  twiddle1 = wavetable->twiddle[i];
	  twiddle2 = twiddle1 + tskip;
	  twiddle3 = twiddle2 + tskip;
	  gsl_fft_halfcomplex_pass_4 (in, istride, out, ostride,
				      product, n, twiddle1, twiddle2, 
				      twiddle3);
	}
      else if (factor == 5)
	{
	  twiddle1 = wavetable->twiddle[i];
	  twiddle2 = twiddle1 + tskip;
	  twiddle3 = twiddle2 + tskip;
	  twiddle4 = twiddle3 + tskip;
	  gsl_fft_halfcomplex_pass_5 (in, istride, out, ostride,
				      product, n, twiddle1, twiddle2, 
				      twiddle3, twiddle4);
	}
      else
	{
	  twiddle1 = wavetable->twiddle[i];
	  gsl_fft_halfcomplex_pass_n (in, istride, out, ostride,
				      factor, product, n, twiddle1);
	}
    }

  if (state == 1)		/* copy results back from scratch to data */
    {
      for (i = 0; i < n; i++)
	{
	  a[istride*i] = b[ostride*i] ;
	}
    }

  return 0;

}
