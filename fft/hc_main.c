#include <config.h>
#include <stdlib.h>
#include <math.h>

#include <gsl_errno.h>
#include <gsl_complex.h>
#include <gsl_fft_halfcomplex.h>

#include "hc_pass.h"

int
gsl_fft_halfcomplex_backward (double data[], const size_t stride, 
			      const size_t n,
			      const gsl_fft_halfcomplex_wavetable * wavetable)
{
  int status = gsl_fft_halfcomplex (data, stride, n, wavetable) ;
  return status ;
}

int
gsl_fft_halfcomplex_inverse (double data[], const size_t stride, 
			     const size_t n,
			     const gsl_fft_halfcomplex_wavetable * wavetable)
{
  int status = gsl_fft_halfcomplex (data, stride, n, wavetable);

  if (status)
    {
      return status;
    }

  /* normalize inverse fft with 1/n */

  {
    const double norm = 1.0 / n;
    size_t i;
    for (i = 0; i < n; i++)
      {
	data[stride*i] *= norm;
      }
  }
  return status;
}

int
gsl_fft_halfcomplex (double data[], const size_t stride, const size_t n,
		     const gsl_fft_halfcomplex_wavetable * wavetable)
{
  size_t factor, product, q, state;
  size_t i;
  size_t nf;
  int product_1;
  int tskip;
  gsl_complex *twiddle1, *twiddle2, *twiddle3, *twiddle4;

  double * const scratch = wavetable->scratch;

  double * in = data;
  size_t istride = stride;
  double * out = scratch;
  size_t ostride = 1;

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
	  in = data;
	  istride = stride;
	  out = scratch;
	  ostride = 1;
	  state = 1;
	}
      else
	{
	  in = scratch;
	  istride = 1;
	  out = data;
	  ostride = stride;
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
	  data[stride*i] = scratch[i] ;
	}
    }

  return 0;

}


