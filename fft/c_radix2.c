#include <config.h>
#include <math.h>

#include <gsl_complex.h>
#include <gsl_fft_complex.h>
#include <gsl_errno.h>

int
gsl_fft_complex_radix2_forward (gsl_complex data[],
				const size_t n)
{
  gsl_fft_direction sign = forward;
  int status = gsl_fft_complex_radix2 (data, n, sign);
  return status;
}

int
gsl_fft_complex_radix2_backward (gsl_complex data[],
				 const size_t n)
{
  gsl_fft_direction sign = backward;
  int status = gsl_fft_complex_radix2 (data, n, sign);
  return status;
}

int
gsl_fft_complex_radix2_inverse (gsl_complex data[],
				const size_t n)
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
    size_t i;
    for (i = 0; i < n; i++)
      {
	data[i].real *= norm;
	data[i].imag *= norm;
      }
  }
  return status;
}


int
gsl_fft_complex_radix2_dif_forward (gsl_complex data[],
				const size_t n)
{
  gsl_fft_direction sign = forward;
  int status = gsl_fft_complex_radix2_dif (data, n, sign);
  return status;
}

int
gsl_fft_complex_radix2_dif_backward (gsl_complex data[],
				 const size_t n)
{
  gsl_fft_direction sign = backward;
  int status = gsl_fft_complex_radix2_dif (data, n, sign);
  return status;
}

int
gsl_fft_complex_radix2_dif_inverse (gsl_complex data[],
				const size_t n)
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
    size_t i;
    for (i = 0; i < n; i++)
      {
	data[i].real *= norm;
	data[i].imag *= norm;
      }
  }
  return status;
}


int
gsl_fft_complex_radix2 (gsl_complex data[],
			const size_t n,
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
  
  status = fft_complex_bitreverse_order(data, n, logn) ;

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
	  
	  const double z1_real = data[j].real ;
	  const double z1_imag = data[j].imag ;

	  const double wd_real = z1_real ;
	  const double wd_imag = z1_imag ;
	  
	  data[j].real = data[i].real - wd_real;
	  data[j].imag = data[i].imag - wd_imag;
	  data[i].real += wd_real;
	  data[i].imag += wd_imag;
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

	      const double z1_real = data[j].real ;
	      const double z1_imag = data[j].imag ;
	      
	      const double wd_real = w_real * z1_real - w_imag * z1_imag;
	      const double wd_imag = w_real * z1_imag + w_imag * z1_real;

	      data[j].real = data[i].real - wd_real;
	      data[j].imag = data[i].imag - wd_imag;
	      data[i].real += wd_real;
	      data[i].imag += wd_imag;
	    }
	}
      dual *= 2;
    }

  return 0;

}



int
gsl_fft_complex_radix2_dif (gsl_complex data[],
			    const size_t n,
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








