#include <config.h>
#include <math.h>
#include <gsl_errno.h>
#include <gsl_fft_halfcomplex.h>

#include "factorize.h"
#include "fft.h"
#include "bitreverse.h"

int
gsl_fft_halfcomplex_radix2_backward (double data[],
				     const size_t stride,
				     const size_t n)
{
  int status = gsl_fft_halfcomplex_radix2 (data, stride, n) ;
  return status ;
}

int
gsl_fft_halfcomplex_radix2_inverse (double data[],
				    const size_t stride,
				    const size_t n)
{
  int status = gsl_fft_halfcomplex_radix2 (data, stride, n);

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
gsl_fft_halfcomplex_radix2 (double data[],
			    const size_t stride,
			    const size_t n)
{
  int result ;
  size_t p, p_1, q;
  size_t i; 
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

  p = n; q = 1 ; p_1 = n/2 ;

  for (i = 1; i <= logn; i++)
    {
      size_t a, b;

      /* a = 0 */

      for (b = 0; b < q; b++)
	{
	  const double z0 = VECTOR(data,stride,b*p);
	  const double z1 = VECTOR(data,stride,b*p + p_1);
	  
	  const double t0_real = z0 + z1 ;
	  const double t1_real = z0 - z1 ;
	  
	  VECTOR(data,stride,b*p) = t0_real;
	  VECTOR(data,stride,b*p + p_1) = t1_real ;
	}

      /* a = 1 ... p_{i-1}/2 - 1 */

      {
	double w_real = 1.0;
	double w_imag = 0.0;

	const double theta = 2.0 * M_PI / p;
	
	const double s = sin (theta);
	const double t = sin (theta / 2.0);
	const double s2 = 2.0 * t * t;
	
	for (a = 1; a < (p_1)/2; a++)
	  {
	    /* trignometric recurrence for w-> exp(i theta) w */
	    
	    {
	      const double tmp_real = w_real - s * w_imag - s2 * w_real;
	      const double tmp_imag = w_imag + s * w_real - s2 * w_imag;
	      w_real = tmp_real;
	      w_imag = tmp_imag;
	    }
	    
	    for (b = 0; b < q; b++)
	      {
		double z0_real = VECTOR(data,stride,b*p + a) ;
		double z0_imag = VECTOR(data,stride,b*p + p - a) ;
		double z1_real = VECTOR(data,stride,b*p + p_1 - a) ;
		double z1_imag = -VECTOR(data,stride,b*p + p_1 + a) ;
		
		/* t0 = z0 + z1 */
		
		double t0_real = z0_real + z1_real;
		double t0_imag = z0_imag + z1_imag;
		
		/* t1 = (z0 - z1) */
		
		double t1_real = z0_real -  z1_real;
		double t1_imag = z0_imag -  z1_imag;
		
		VECTOR(data,stride,b*p + a) = t0_real ;
		VECTOR(data,stride,b*p + p_1 - a) = t0_imag ;
		
		VECTOR(data,stride,b*p + p_1 + a) = (w_real * t1_real - w_imag * t1_imag) ;
		VECTOR(data,stride,b*p + p - a) = (w_real * t1_imag + w_imag * t1_real) ;
	      }
	  }
      }

      if (p_1 >  1) {
	for (b = 0; b < q; b++) {
	  VECTOR(data,stride,b*p + p_1/2) *= 2 ;
	  VECTOR(data,stride,b*p + p_1 + p_1/2) *= -2 ;
	}
      }

      p_1 = p_1 / 2 ;
      p = p / 2 ;
      q = q * 2 ;
    }

  /* bit reverse the ordering of output data for decimation in
     frequency algorithm */
  
  status = fft_real_bitreverse_order(data, stride, n, logn) ;

  return 0;

}
