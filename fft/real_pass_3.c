#include <config.h>
#include <math.h>

#include <gsl_complex.h>
#include <gsl_fft_real.h>

#include "real_pass.h"

int
gsl_fft_real_pass_3 (const double in[], 
		     const size_t istride,
		     double out[],
		     const size_t ostride,
		     const size_t product,
		     const size_t n,
		     const gsl_complex twiddle1[],
		     const gsl_complex twiddle2[])
{
  size_t k, k1;

  const size_t factor = 3;
  const size_t m = n / factor;
  const size_t q = n / product;
  const size_t product_1 = product / factor;

  const double tau = sqrt (3.0) / 2.0;

  for (k1 = 0; k1 < q; k1++)
    {
      const size_t from0 = k1 * product_1;
      const size_t from1 = from0 + m;
      const size_t from2 = from1 + m;
      
      const double z0_real = VECTOR(in,istride,from0);
      const double z1_real = VECTOR(in,istride,from1);
      const double z2_real = VECTOR(in,istride,from2);
      
      const double t1 = z1_real + z2_real;

      const double x0_real = z0_real + t1;
      const double x1_real = z0_real - t1 / 2.0;
      const double x1_imag = -tau * (z1_real - z2_real);

      const size_t to0 = product * k1;
      const size_t to1 = to0 + 2 * product_1 - 1;

      VECTOR(out,ostride,to0) = x0_real;
      VECTOR(out,ostride,to1) = x1_real;
      VECTOR(out,ostride,to1 + 1) = x1_imag;
    }

  if (product_1 == 1)
    return 0;

  for (k = 1; k < (product_1 + 1) / 2; k++)
    {
      const double w1_real = GSL_REAL(twiddle1[k - 1]);
      const double w1_imag = -GSL_IMAG(twiddle1[k - 1]);
      const double w2_real = GSL_REAL(twiddle2[k - 1]);
      const double w2_imag = -GSL_IMAG(twiddle2[k - 1]);

      for (k1 = 0; k1 < q; k1++)
	{
	  const size_t from0 = k1 * product_1 + 2 * k - 1;
	  const size_t from1 = from0 + m;
	  const size_t from2 = from1 + m;
	  
	  const double f0_real = VECTOR(in,istride,from0);
	  const double f0_imag = VECTOR(in,istride,from0 + 1);
	  const double f1_real = VECTOR(in,istride,from1);
	  const double f1_imag = VECTOR(in,istride,from1 + 1);
	  const double f2_real = VECTOR(in,istride,from2);
	  const double f2_imag = VECTOR(in,istride,from2 + 1);
	  
	  const double z0_real = f0_real;
	  const double z0_imag = f0_imag;
	  const double z1_real = w1_real * f1_real - w1_imag * f1_imag;
	  const double z1_imag = w1_real * f1_imag + w1_imag * f1_real;
	  const double z2_real = w2_real * f2_real - w2_imag * f2_imag;
	  const double z2_imag = w2_real * f2_imag + w2_imag * f2_real;

	  /* compute x = W(3) z */

	  /* t1 = z1 + z2 */
	  const double t1_real = z1_real + z2_real;
	  const double t1_imag = z1_imag + z2_imag;
	  
	  /* t2 = z0 - t1/2 */
	  const double t2_real = z0_real - t1_real / 2;
	  const double t2_imag = z0_imag - t1_imag / 2;
	  
	  /* t3 = (+/-) sin(pi/3)*(z1 - z2) */
	  const double t3_real = -tau * (z1_real - z2_real);
	  const double t3_imag = -tau * (z1_imag - z2_imag);
	  
	  /* x0 = z0 + t1 */
	  const double x0_real = z0_real + t1_real;
	  const double x0_imag = z0_imag + t1_imag;
	  
	  /* x1 = t2 + i t3 */
	  const double x1_real = t2_real - t3_imag;
	  const double x1_imag = t2_imag + t3_real;

	    /* x2 = t2 - i t3 */
	  const double x2_real = t2_real + t3_imag;
	  const double x2_imag = t2_imag - t3_real;

	  /* apply twiddle factors */
	  
	  const size_t to0 = k1 * product + 2 * k - 1;
	  const size_t to1 = to0 + 2 * product_1;
	  const size_t to2 = 2 * product_1 - 2 * k + k1 * product - 1;
	  
	  /* to0 = 1 * x0 */
	  VECTOR(out,ostride,to0) = x0_real;
	  VECTOR(out,ostride,to0 + 1) = x0_imag;
	  
	  /* to1 = 1 * x1 */
	  VECTOR(out,ostride,to1) = x1_real;
	  VECTOR(out,ostride,to1 + 1) = x1_imag;
	  
	  /* to2 = 1 * x2 */
	  VECTOR(out,ostride,to2) = x2_real;
	  VECTOR(out,ostride,to2 + 1) = -x2_imag;
	}
    }
  
  if (product_1 % 2 == 1)
    return 0;

  for (k1 = 0; k1 < q; k1++)
    {
      const size_t from0 = k1 * product_1 + product_1 - 1;
      const size_t from1 = from0 + m;
      const size_t from2 = from1 + m;
      
      const double z0_real = VECTOR(in,istride,from0);
      const double z1_real = VECTOR(in,istride,from1);
      const double z2_real = VECTOR(in,istride,from2);

      const double t1 = z1_real - z2_real;
      const double x0_real = z0_real + t1 / 2.0;
      const double x0_imag = -tau * (z1_real + z2_real);
      const double x1_real = z0_real - t1;

      const size_t to0 = k1 * product + product_1 - 1;
      const size_t to1 = to0 + 2 * product_1;
      
      VECTOR(out,ostride,to0) = x0_real;
      VECTOR(out,ostride,to0 + 1) = x0_imag;
      VECTOR(out,ostride,to1) = x1_real;
    }

  return 0;
}
