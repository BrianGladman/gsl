#include <config.h>
#include <stdlib.h>
#include <math.h>

#include <gsl_errno.h>
#include <gsl_complex.h>
#include <gsl_fft_complex.h>

#include "fft_complex.h"

int
gsl_fft_complex_forward (gsl_vector_complex * data,
			 const gsl_fft_complex_wavetable * wavetable)
{
  gsl_fft_direction sign = forward;
  int status = gsl_fft_complex (data, wavetable, sign);
  return status;
}

int
gsl_fft_complex_backward (gsl_vector_complex * data,
			  const gsl_fft_complex_wavetable * wavetable)
{
  gsl_fft_direction sign = backward;
  int status = gsl_fft_complex (data, wavetable, sign);
  return status;
}

int
gsl_fft_complex_inverse (gsl_vector_complex * data,
			 const gsl_fft_complex_wavetable * wavetable)
{
  gsl_fft_direction sign = backward;
  int status = gsl_fft_complex (data, wavetable, sign);

  if (status)
    {
      return status;
    }

  const size_t n = data->size ;
  const size_t stride = data->stride ;
  double * const v = data->data ;

  /* normalize inverse fft with 1/n */

  {
    const double norm = 1.0 / n;
    size_t i;
    for (i = 0; i < n; i++)
      {
	REAL(data,stride,i) *= norm;
	IMAG(data,stride,i) *= norm;
      }
  }
  return status;
}

int
gsl_fft_complex (gsl_vector_complex * data,
		 const gsl_fft_complex_wavetable * wavetable,
		 const gsl_fft_direction sign)
{

  const size_t nf = wavetable->nf;

  size_t i;

  size_t q, product = 1;

  gsl_complex *scratch = wavetable->scratch;
  gsl_complex *twiddle1, *twiddle2, *twiddle3, *twiddle4, *twiddle5, *twiddle6;

  size_t state = 0;

  const size_t n = data->size ;

  double * const a = data->data;
  double * const b = scratch;

  const size_t astride = data->stride ;
  const size_t bstride = 1 ;

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
      const size_t factor = wavetable->factor[i];
      product *= factor;
      q = n / product;

      if (state == 0)
	{
	  in = a;
	  istride = astride;
	  out = b;
	  ostrid = bstride;
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
	  gsl_fft_complex_pass_2 (in, istride, out, ostride, sign, 
				  product, n, twiddle1);
	}
      else if (factor == 3)
	{
	  twiddle1 = wavetable->twiddle[i];
	  twiddle2 = twiddle1 + q;
	  gsl_fft_complex_pass_3 (in, istride, out, ostride, sign, 
				  product, n, twiddle1, twiddle2);
	}
      else if (factor == 4)
	{
	  twiddle1 = wavetable->twiddle[i];
	  twiddle2 = twiddle1 + q;
	  twiddle3 = twiddle2 + q;
	  gsl_fft_complex_pass_4 (in, istride, out, ostride, sign, 
				  product, n, twiddle1, twiddle2, twiddle3);
	}
      else if (factor == 5)
	{
	  twiddle1 = wavetable->twiddle[i];
	  twiddle2 = twiddle1 + q;
	  twiddle3 = twiddle2 + q;
	  twiddle4 = twiddle3 + q;
	  gsl_fft_complex_pass_5 (in, istride, out, ostride, sign, 
				  product, n, twiddle1, twiddle2, twiddle3, 
				  twiddle4);
	}
      else if (factor == 6)
	{
	  twiddle1 = wavetable->twiddle[i];
	  twiddle2 = twiddle1 + q;
	  twiddle3 = twiddle2 + q;
	  twiddle4 = twiddle3 + q;
	  twiddle5 = twiddle4 + q;
	  gsl_fft_complex_pass_6 (in, istride, out, ostride, sign, 
				  product, n, twiddle1, twiddle2, twiddle3, 
				  twiddle4, twiddle5);
	}
      else if (factor == 7)
	{
	  twiddle1 = wavetable->twiddle[i];
	  twiddle2 = twiddle1 + q;
	  twiddle3 = twiddle2 + q;
	  twiddle4 = twiddle3 + q;
	  twiddle5 = twiddle4 + q;
	  twiddle6 = twiddle5 + q;
	  gsl_fft_complex_pass_7 (in, istride, out, ostride, sign, 
				  product, n, twiddle1, twiddle2, twiddle3, 
				  twiddle4, twiddle5, twiddle6);
	}
      else
	{
	  twiddle1 = wavetable->twiddle[i];
	  gsl_fft_complex_pass_n (in, istride, out, ostride, sign, 
				  factor, product, n, twiddle1);
	}
    }

  if (state == 1)		/* copy results back from scratch to data */
    {
      for (i = 0; i < n; i++)
	{
	  REAL(a,istride,i) = REAL(b,ostride,i) ;
	  IMAG(a,istride,i) = IMAG(b,ostride,i) ;
	}
    }

  return 0;

}
