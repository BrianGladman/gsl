#include <config.h>
#include <stdlib.h>
#include <math.h>

#include <gsl_errno.h>
#include <gsl_complex.h>
#include <gsl_fft_real.h>

#include "fft_real.h"

int
gsl_fft_real (double data[], const size_t stride, const size_t n,
	      const gsl_fft_real_wavetable * wavetable)
{
  const size_t nf = wavetable->nf;

  size_t i;

  size_t q, product = 1;
  size_t tskip;
  size_t product_1;

  double *scratch = wavetable->scratch;
  gsl_complex *twiddle1, *twiddle2, *twiddle3, *twiddle4;

  size_t state = 0;
  double *in = data;
  size_t istride = stride ;
  double *out = scratch;
  size_t ostride = 1 ;
  
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

  for (i = 0; i < nf; i++)
    {
      const size_t factor = wavetable->factor[i];
      product_1 = product;
      product *= factor;
      q = n / product;

      tskip = (product_1 + 1) / 2 - 1;

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
	  gsl_fft_real_pass_2 (in, istride, out, ostride, product, n, twiddle1);
	}
      else if (factor == 3)
	{
	  twiddle1 = wavetable->twiddle[i];
	  twiddle2 = twiddle1 + tskip;
	  gsl_fft_real_pass_3 (in, istride, out, ostride, product, n, twiddle1,
			       twiddle2);
	}
      else if (factor == 4)
	{
	  twiddle1 = wavetable->twiddle[i];
	  twiddle2 = twiddle1 + tskip;
	  twiddle3 = twiddle2 + tskip;
	  gsl_fft_real_pass_4 (in, istride, out, ostride, product, n, twiddle1,
			       twiddle2, twiddle3);
	}
      else if (factor == 5)
	{
	  twiddle1 = wavetable->twiddle[i];
	  twiddle2 = twiddle1 + tskip;
	  twiddle3 = twiddle2 + tskip;
	  twiddle4 = twiddle3 + tskip;
	  gsl_fft_real_pass_5 (in, istride, out, ostride, product, n, twiddle1,
			       twiddle2, twiddle3, twiddle4);
	}
      else
	{
	  twiddle1 = wavetable->twiddle[i];
	  gsl_fft_real_pass_n (in, istride, out, ostride, factor, product, n,
			       twiddle1);
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
