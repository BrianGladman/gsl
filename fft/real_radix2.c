#include "complex_internal.h"
#include "factorize.h"
#include "bitreverse.h"

int
FUNCTION(gsl_fft_real,radix2_transform) (BASE data[], const size_t stride,  const size_t n)
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

  result = fft_binary_logn(n) ;

  if (result == -1) 
    {
      GSL_ERROR ("n is not a power of 2", GSL_EINVAL);
    } 
  else 
    {
      logn = result ;
    }

  /* bit reverse the ordering of input data for decimation in time algorithm */
  
  status = FUNCTION(fft_real,bitreverse_order)(data, stride, n, logn) ;

  /* apply fft recursion */

  p = 1; q = n ;

  for (i = 1; i <= logn; i++)
    {
      size_t a, b;

      p_1 = p ;
      p = 2 * p ;
      q = q / 2 ;

      /* a = 0 */

      for (b = 0; b < q; b++)
	{
	  double t0_real = VECTOR(data,stride,b*p) + VECTOR(data,stride,b*p + p_1) ;
	  double t1_real = VECTOR(data,stride,b*p) - VECTOR(data,stride,b*p + p_1) ;
	  
	  VECTOR(data,stride,b*p) = t0_real ;
	  VECTOR(data,stride,b*p + p_1) = t1_real ;
	}

      /* a = 1 ... p_{i-1}/2 - 1 */

      {
	double w_real = 1.0;
	double w_imag = 0.0;

	const double theta = - 2.0 * M_PI / p;
	
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
		double z0_imag = VECTOR(data,stride,b*p + p_1 - a) ;
		double z1_real = VECTOR(data,stride,b*p + p_1 + a) ;
		double z1_imag = VECTOR(data,stride,b*p + p - a) ;
		
		/* t0 = z0 + w * z1 */
		
		double t0_real = z0_real + w_real * z1_real - w_imag * z1_imag;
		double t0_imag = z0_imag + w_real * z1_imag + w_imag * z1_real;
		
		/* t1 = z0 - w * z1 */
		
		double t1_real = z0_real - w_real * z1_real + w_imag * z1_imag;
		double t1_imag = z0_imag - w_real * z1_imag - w_imag * z1_real;
		
		VECTOR(data,stride,b*p + a) = t0_real ;
		VECTOR(data,stride,b*p + p - a) = t0_imag ;
		
		VECTOR(data,stride,b*p + p_1 - a) = t1_real ;
		VECTOR(data,stride,b*p + p_1 + a) = -t1_imag ;
	      }
	  }
      }

      if (p_1 >  1) 
	{
	  for (b = 0; b < q; b++) 
	    {
	      /* a = p_{i-1}/2 */
	      
	      VECTOR(data,stride,b*p + p - p_1/2) *= -1 ;
	    }
	}
    }
  return 0;
}
