#include <config.h>
#include <math.h>

#include <gsl_complex.h>
#include <gsl_fft_complex.h>
#include <gsl_errno.h>

#include "fft.h"
#include "bitreverse.h"

int
gsl_fft_complex_radix2_forward (double data[], const size_t stride, const size_t n)
{
  gsl_fft_direction sign = forward;
  int status = gsl_fft_complex_radix2 (data, stride, n, sign);
  return status;
}

int
gsl_fft_complex_radix2_backward (double data[], const size_t stride, const size_t n)
{
  gsl_fft_direction sign = backward;
  int status = gsl_fft_complex_radix2 (data, stride, n, sign);
  return status;
}

int
gsl_fft_complex_radix2_inverse (double data[], const size_t stride, const size_t n)
{
  gsl_fft_direction sign = backward;
  int status = gsl_fft_complex_radix2 (data, stride, n, sign);

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
gsl_fft_complex_radix2_dif_forward (double data[], const size_t stride, const size_t n)
{
  gsl_fft_direction sign = forward;
  int status = gsl_fft_complex_radix2_dif (data, stride, n, sign);
  return status;
}

int
gsl_fft_complex_radix2_dif_backward (double data[], const size_t stride, const size_t n)
{
  gsl_fft_direction sign = backward;
  int status = gsl_fft_complex_radix2_dif (data, stride, n, sign);
  return status;
}

int
gsl_fft_complex_radix2_dif_inverse (double data[], const size_t stride, const size_t n)
{
  gsl_fft_direction sign = backward;
  int status = gsl_fft_complex_radix2_dif (data, stride, n, sign);

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
gsl_fft_complex_radix2 (double data[], const size_t stride, const size_t n,
			const gsl_fft_direction sign)
{

  int result ;
  size_t dual;
  size_t bit; 
  size_t logn = 0;
  int status;

  if (n == 1) /* identity operation */
    {
      return 0 ;
    }

  /* make sure that n is a power of 2 */

  result = gsl_fft_binary_logn(n) ;

  if (result == -1) 
    {
      GSL_ERROR ("n is not a power of 2", GSL_EINVAL);
    } 
  else 
    {
      logn = result ;
    }

  /* bit reverse the ordering of input data for decimation in time algorithm */
  
  status = fft_complex_bitreverse_order(data, stride, n, logn) ;

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

      size_t a, b;

      /* a = 0 */

      for (b = 0; b < n; b += 2 * dual)
	{
	  const size_t i = b ;
	  const size_t j = b + dual;
	  
	  const double z1_real = REAL(data,stride,j) ;
	  const double z1_imag = IMAG(data,stride,j) ;

	  const double wd_real = z1_real ;
	  const double wd_imag = z1_imag ;
	  
	  REAL(data,stride,j) = REAL(data,stride,i) - wd_real;
	  IMAG(data,stride,j) = IMAG(data,stride,i) - wd_imag;
	  REAL(data,stride,i) += wd_real;
	  IMAG(data,stride,i) += wd_imag;
	}
      
      /* a = 1 .. (dual-1) */

      for (a = 1; a < dual; a++)
	{

	  /* trignometric recurrence for w-> exp(i theta) w */

	  {
	    const double tmp_real = w_real - s * w_imag - s2 * w_real;
	    const double tmp_imag = w_imag + s * w_real - s2 * w_imag;
	    w_real = tmp_real;
	    w_imag = tmp_imag;
	  }

	  for (b = 0; b < n; b += 2 * dual)
	    {
	      const size_t i = b + a;
	      const size_t j = b + a + dual;

	      const double z1_real = REAL(data,stride,j) ;
	      const double z1_imag = IMAG(data,stride,j) ;
	      
	      const double wd_real = w_real * z1_real - w_imag * z1_imag;
	      const double wd_imag = w_real * z1_imag + w_imag * z1_real;

	      REAL(data,stride,j) = REAL(data,stride,i) - wd_real;
	      IMAG(data,stride,j) = IMAG(data,stride,i) - wd_imag;
	      REAL(data,stride,i) += wd_real;
	      IMAG(data,stride,i) += wd_imag;
	    }
	}
      dual *= 2;
    }

  return 0;

}



int
gsl_fft_complex_radix2_dif (double data[], const size_t stride, const size_t n,
			    const gsl_fft_direction sign)
{
  int result ;
  size_t dual;
  size_t bit; 
  size_t logn = 0;
  int status;

  if (n == 1) /* identity operation */
    {
      return 0 ;
    }

  /* make sure that n is a power of 2 */

  result = gsl_fft_binary_logn(n) ;

  if (result == -1) 
    {
      GSL_ERROR ("n is not a power of 2", GSL_EINVAL);
    } 
  else 
    {
      logn = result ;
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

      size_t a, b;

      for (b = 0; b < dual; b++)
	{
	  for (a = 0; a < n; a+= 2 * dual)
	    {
	      const size_t i = b + a;
	      const size_t j = b + a + dual;
	      
	      const double t1_real = REAL(data,stride,i) + REAL(data,stride,j);
	      const double t1_imag = IMAG(data,stride,i) + IMAG(data,stride,j);
	      const double t2_real = REAL(data,stride,i) - REAL(data,stride,j);
	      const double t2_imag = IMAG(data,stride,i) - IMAG(data,stride,j);

	      REAL(data,stride,i) = t1_real;
	      IMAG(data,stride,i) = t1_imag;
	      REAL(data,stride,j) = w_real*t2_real - w_imag * t2_imag;
	      IMAG(data,stride,j) = w_real*t2_imag + w_imag * t2_real;
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
  
  status = fft_complex_bitreverse_order(data, stride, n, logn) ;

  return 0;

}








