#include <math.h>
#include <gsl_complex.h>
#include <gsl_fft_complex.h>
#include <gsl_errno.h>

int
gsl_fft_complex_radix2_forward (complex data[],
				const unsigned int n)
{
  gsl_fft_direction sign = forward;
  int status = gsl_fft_complex_radix2 (data, n, sign);
  return status;
}

int
gsl_fft_complex_radix2_backward (complex data[],
				 const unsigned int n)
{
  gsl_fft_direction sign = backward;
  int status = gsl_fft_complex_radix2 (data, n, sign);
  return status;
}

int
gsl_fft_complex_radix2_inverse (complex data[],
				const unsigned int n)
{
  gsl_fft_direction sign = backward;
  int status = gsl_fft_complex_radix2 (data, n, sign);

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
gsl_fft_complex_radix2_dif_forward (complex data[],
				const unsigned int n)
{
  gsl_fft_direction sign = forward;
  int status = gsl_fft_complex_radix2_dif (data, n, sign);
  return status;
}

int
gsl_fft_complex_radix2_dif_backward (complex data[],
				 const unsigned int n)
{
  gsl_fft_direction sign = backward;
  int status = gsl_fft_complex_radix2_dif (data, n, sign);
  return status;
}

int
gsl_fft_complex_radix2_dif_inverse (complex data[],
				const unsigned int n)
{
  gsl_fft_direction sign = backward;
  int status = gsl_fft_complex_radix2_dif (data, n, sign);

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
gsl_fft_complex_radix2 (complex data[],
			const unsigned int n,
			const gsl_fft_direction sign)
{

  int dual;
  int bit; unsigned int logn = 0;
  int status;

  if (n == 1) /* identity operation */
    {
      return 0 ;
    }

  /* make sure that n is a power of 2 */

  logn = gsl_fft_binary_logn(n) ;

  if (logn == -1) {
    GSL_ERROR ("n is not a power of 2", GSL_EINVAL);
    return -1;
  }

  /* bit reverse the ordering of input data for decimation in time algorithm */
  
  status = gsl_fft_complex_bitreverse_order(data, n, logn) ;

  /* apply fft recursion */

  dual = 1;

  for (bit = 0; bit < logn; bit++)
    {
      double w_real = 1.0;
      double w_imag = 0.0;

      const double theta = 2.0 * ((int) sign) * M_PI / (2.0 * (double) dual);

      const double s = sin (theta);
      const double t = sin (theta / 2.0);
      const double s2 = 2.0 * t * t;

      unsigned int a, b;

      for (a = 0; a < dual; a++)
	{
	  for (b = 0; b < n; b += 2 * dual)
	    {
	      const unsigned int i = b + a;
	      const unsigned int j = b + a + dual;

	      const double wd_real = w_real * data[j].real - w_imag * data[j].imag;
	      const double wd_imag = w_real * data[j].imag + w_imag * data[j].real;

	      data[j].real = data[i].real - wd_real;
	      data[j].imag = data[i].imag - wd_imag;
	      data[i].real += wd_real;
	      data[i].imag += wd_imag;
	    }

	  /* trignometric recurrence for w-> exp(i theta) w */

	  {
	    const double tmp_real = w_real - s * w_imag - s2 * w_real;
	    const double tmp_imag = w_imag + s * w_real - s2 * w_imag;
	    w_real = tmp_real;
	    w_imag = tmp_imag;
	  }
	}
      dual *= 2;
    }

  return 0;

}



int
gsl_fft_complex_radix2_dif (complex data[],
			    const unsigned int n,
			    const gsl_fft_direction sign)
{
  int dual;
  int bit; 
  unsigned int logn = 0;
  int status;

  if (n == 1) /* identity operation */
    {
      return 0 ;
    }

  /* make sure that n is a power of 2 */

  logn = gsl_fft_binary_logn(n) ;

  if (logn == -1) {
    GSL_ERROR ("n is not a power of 2", GSL_EINVAL);
    return -1;
  }

  /* apply fft recursion */

  dual = n / 2;

  for (bit = 0; bit < logn; bit++)
    {
      double w_real = 1.0;
      double w_imag = 0.0;

      const double theta = 2.0 * ((int) sign) * M_PI / ((double) (2 * dual));

      const double s = sin (theta);
      const double t = sin (theta / 2.0);
      const double s2 = 2.0 * t * t;

      unsigned int a, b;

      for (b = 0; b < dual; b++)
	{
	  for (a = 0; a < n; a+= 2 * dual)
	    {
	      const unsigned int i = b + a;
	      const unsigned int j = b + a + dual;
	      
	      const double t1_real = data[i].real + data[j].real;
	      const double t1_imag = data[i].imag + data[j].imag;
	      const double t2_real = data[i].real - data[j].real;
	      const double t2_imag = data[i].imag - data[j].imag;

	      data[i].real = t1_real;
	      data[i].imag = t1_imag;
	      data[j].real = w_real*t2_real - w_imag * t2_imag;
	      data[j].imag = w_real*t2_imag + w_imag * t2_real;
	    }

	  /* trignometric recurrence for w-> exp(i theta) w */

	  {
	    const double tmp_real = w_real - s * w_imag - s2 * w_real;
	    const double tmp_imag = w_imag + s * w_real - s2 * w_imag;
	    w_real = tmp_real;
	    w_imag = tmp_imag;
	  }
	}
      dual /= 2;
    }

  /* bit reverse the ordering of output data for decimation in
     frequency algorithm */
  
  status = gsl_fft_complex_bitreverse_order(data, n, logn) ;

  return 0;

}


int gsl_fft_binary_logn (const unsigned int n)
{
  int binary_logn = 0 ;
  unsigned int k = 1;

  while (k < n)
    {
      k *= 2;
      binary_logn++;
    }

  if (n != (1 << binary_logn))
    {
      /* n is not a power of 2 */
      return -1 ; 
    } 
  else 
    {
      return binary_logn;
    }
      
}


int gsl_fft_complex_bitreverse_order (complex data[], 
				      const unsigned int n,
				      const unsigned int logn)
{
  unsigned int i;

  for (i = 0; i < n; i++)
    {
      unsigned int j = 0;
      unsigned int i_tmp = i;
      unsigned int bit;

      for (bit = 0; bit < logn; bit++)
	{
	  j <<= 1;		/* reverse shift i into j */
	  j |= i_tmp & 1;
	  i_tmp >>= 1;
	}
      
      if (i < j)
	{
	  const complex data_tmp = data[i];
	  data[i] = data[j];
	  data[j] = data_tmp;
	}
    }
  return 0;
}



