int
FUNCTION(gsl_fft_complex,forward) (TYPE(gsl_complex_packed_array) data, 
				   const size_t stride, 
				   const size_t n,
				   const TYPE(gsl_fft_wavetable_complex) * wavetable)
{
  gsl_fft_direction sign = forward;
  int status = FUNCTION(gsl_fft_complex,transform) (data, stride, n, 
						    wavetable, sign);
  return status;
}

int
FUNCTION(gsl_fft_complex,backward) (TYPE(gsl_complex_packed_array) data,
				    const size_t stride, 
				    const size_t n,
				    const TYPE(gsl_fft_wavetable_complex) * wavetable)
{
  gsl_fft_direction sign = backward;
  int status = FUNCTION(gsl_fft_complex,transform) (data, stride, n, 
						    wavetable, sign);
  return status;
}

int
FUNCTION(gsl_fft_complex,inverse) (TYPE(gsl_complex_packed_array) data, 
				   const size_t stride, 
				   const size_t n,
				   const TYPE(gsl_fft_wavetable_complex) * wavetable)
{
  gsl_fft_direction sign = backward;
  int status = FUNCTION(gsl_fft_complex,transform) (data, stride, n, 
						    wavetable, sign);

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
	REAL(data,stride,i) *= norm;
	IMAG(data,stride,i) *= norm;
      }
  }
  return status;
}

int
FUNCTION(gsl_fft_complex,transform) (TYPE(gsl_complex_packed_array) data, 
				     const size_t stride, 
				     const size_t n,
				     const TYPE(gsl_fft_wavetable_complex) * wavetable,
				     const gsl_fft_direction sign)
{
  const size_t nf = wavetable->nf;

  size_t i;

  size_t q, product = 1;

  gsl_complex *twiddle1, *twiddle2, *twiddle3, *twiddle4,
    *twiddle5, *twiddle6;

  size_t state = 0;

  BASE * const scratch = wavetable->scratch;

  BASE * in = data;
  size_t istride = stride;

  BASE * out = scratch;
  size_t ostride = 1;

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
	  FUNCTION(fft_complex,pass_2) (in, istride, out, ostride, sign, 
					product, n, twiddle1);
	}
      else if (factor == 3)
	{
	  twiddle1 = wavetable->twiddle[i];
	  twiddle2 = twiddle1 + q;
	  FUNCTION(fft_complex,pass_3) (in, istride, out, ostride, sign, 
					product, n, twiddle1, twiddle2);
	}
      else if (factor == 4)
	{
	  twiddle1 = wavetable->twiddle[i];
	  twiddle2 = twiddle1 + q;
	  twiddle3 = twiddle2 + q;
	  FUNCTION(fft_complex,pass_4) (in, istride, out, ostride, sign, 
					product, n, twiddle1, twiddle2, 
					twiddle3);
	}
      else if (factor == 5)
	{
	  twiddle1 = wavetable->twiddle[i];
	  twiddle2 = twiddle1 + q;
	  twiddle3 = twiddle2 + q;
	  twiddle4 = twiddle3 + q;
	  FUNCTION(fft_complex,pass_5) (in, istride, out, ostride, sign, 
					product, n, twiddle1, twiddle2, 
					twiddle3, twiddle4);
	}
      else if (factor == 6)
	{
	  twiddle1 = wavetable->twiddle[i];
	  twiddle2 = twiddle1 + q;
	  twiddle3 = twiddle2 + q;
	  twiddle4 = twiddle3 + q;
	  twiddle5 = twiddle4 + q;
	  FUNCTION(fft_complex,pass_6) (in, istride, out, ostride, sign, 
					product, n, twiddle1, twiddle2, 
					twiddle3, twiddle4, twiddle5);
	}
      else if (factor == 7)
	{
	  twiddle1 = wavetable->twiddle[i];
	  twiddle2 = twiddle1 + q;
	  twiddle3 = twiddle2 + q;
	  twiddle4 = twiddle3 + q;
	  twiddle5 = twiddle4 + q;
	  twiddle6 = twiddle5 + q;
	  FUNCTION(fft_complex,pass_7) (in, istride, out, ostride, sign, 
					product, n, twiddle1, twiddle2, 
					twiddle3, twiddle4, twiddle5, 
					twiddle6);
	}
      else
	{
	  twiddle1 = wavetable->twiddle[i];
	  FUNCTION(fft_complex,pass_n) (in, istride, out, ostride, sign, 
					factor, product, n, twiddle1);
	}
    }

  if (state == 1)		/* copy results back from scratch to data */
    {
      for (i = 0; i < n; i++)
	{
	  REAL(data,stride,i) = REAL(scratch,1,i) ;
	  IMAG(data,stride,i) = IMAG(scratch,1,i) ;
	}
    }

  return 0;

}
