#include "real_pass.h"
 
void
FUNCTION(fft_real,pass_2) (const BASE in[],
			   const size_t istride,
			   BASE out[],
			   const size_t ostride,
			   const size_t product,
			   const size_t n,
			   const gsl_complex twiddle[])
{
  size_t k, k1;

  const size_t factor = 2;
  const size_t m = n / factor;
  const size_t q = n / product;
  const size_t product_1 = product / factor;

  for (k1 = 0; k1 < q; k1++)
    {
      const size_t from0 = k1 * product_1;
      const size_t from1 = from0 + m;

      const double r0 = VECTOR(in,istride,from0);
      const double r1 = VECTOR(in,istride,from1);
      
      const double s0 = r0 + r1;
      const double s1 = r0 - r1;
      
      const size_t to0 = product * k1;
      const size_t to1 = to0 + product - 1;
      
      VECTOR(out,ostride,to0) = s0;
      VECTOR(out,ostride,to1) = s1;
    }

  if (product_1 == 1)
    return;

  for (k = 1; k < (product_1 + 1) / 2; k++)
    {

      /* forward real transform: w -> conjugate(w) */
      const double w_real = GSL_REAL(twiddle[k - 1]);
      const double w_imag = -GSL_IMAG(twiddle[k - 1]);

      for (k1 = 0; k1 < q; k1++)
	{
	  const size_t from0 = k1 * product_1 + 2 * k - 1;
	  const size_t from1 = from0 + m;

	  const double f0_real = VECTOR(in,istride,from0);
	  const double f0_imag = VECTOR(in,istride,from0 + 1);

	  const double f1_real = VECTOR(in,istride,from1);
	  const double f1_imag = VECTOR(in,istride,from1 + 1);

	  const double z0_real = f0_real;
	  const double z0_imag = f0_imag;

	  const double z1_real = w_real * f1_real - w_imag * f1_imag;
	  const double z1_imag = w_real * f1_imag + w_imag * f1_real;

	  /* compute x = W(2) z */

	  /* x0 = z0 + z1 */
	  const double x0_real = z0_real + z1_real;
	  const double x0_imag = z0_imag + z1_imag;

	  /* x1 = z0 - z1 */
	  const double x1_real = z0_real - z1_real;
	  const double x1_imag = z0_imag - z1_imag;

	  const size_t to0 = k1 * product + 2 * k - 1;
	  const size_t to1 = k1 * product + product - 2 * k - 1;
	  
	  VECTOR(out,ostride,to0) = x0_real;
	  VECTOR(out,ostride,to0 + 1) = x0_imag;
	  
	  /* stored in conjugate location */
	  VECTOR(out,ostride,to1) = x1_real;
	  VECTOR(out,ostride,to1 + 1) = -x1_imag;
	}
    }
  
  if (product_1 % 2 == 1)
    return;

  for (k1 = 0; k1 < q; k1++)
    {
      const size_t from0 = k1 * product_1 + product_1 - 1;
      const size_t from1 = from0 + m;
      const size_t to0 = k1 * product + product_1 - 1;

      VECTOR(out,ostride,to0) = VECTOR(in,istride,from0);
      VECTOR(out,ostride,to0 + 1) = -VECTOR(in,istride,from1);
    }
  return;
}
