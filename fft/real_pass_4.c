#include "real_pass.h"
 
void
FUNCTION(fft_real,pass_4) (const BASE in[],
			   const size_t istride,
			   BASE out[],
			   const size_t ostride,
			   const size_t product,
			   const size_t n,
			   const gsl_complex twiddle1[],
			   const gsl_complex twiddle2[],
			   const gsl_complex twiddle3[])
{
  size_t k, k1;

  const size_t factor = 4;
  const size_t m = n / factor;
  const size_t q = n / product;
  const size_t product_1 = product / factor;

  for (k1 = 0; k1 < q; k1++)
    {
      const size_t from0 = k1 * product_1;
      const size_t from1 = from0 + m;
      const size_t from2 = from1 + m;
      const size_t from3 = from2 + m;
      
      const double z0_real = VECTOR(in,istride,from0);
      const double z1_real = VECTOR(in,istride,from1);
      const double z2_real = VECTOR(in,istride,from2);
      const double z3_real = VECTOR(in,istride,from3);

      /* compute x = W(4) z */

      /* t1 = z0 + z2 */
      const double t1_real = z0_real + z2_real;
      
      /* t2 = z1 + z3 */
      const double t2_real = z1_real + z3_real;
      
	/* t3 = z0 - z2 */
      const double t3_real = z0_real - z2_real;
      
      /* t4 = - (z1 - z3) */
      const double t4_real = -(z1_real - z3_real);
      
      /* x0 = t1 + t2 */
      const double x0_real = t1_real + t2_real;

      /* x1 = t3 + i t4 */
      const double x1_real = t3_real;
      const double x1_imag = t4_real;

      /* x2 = t1 - t2 */
      const double x2_real = t1_real - t2_real;

      const size_t to0 = product * k1;
      const size_t to1 = to0 + 2 * product_1 - 1;
      const size_t to2 = to1 + 2 * product_1;
      
      VECTOR(out,ostride,to0) = x0_real;
      VECTOR(out,ostride,to1) = x1_real;
      VECTOR(out,ostride,to1 + 1) = x1_imag;
      VECTOR(out,ostride,to2) = x2_real;
    }

  if (product_1 == 1)
    return;

  for (k = 1; k < (product_1 + 1) / 2; k++)
    {
      double w1_real, w1_imag, w2_real, w2_imag, w3_real, w3_imag;
      w1_real = GSL_REAL(twiddle1[k - 1]);
      w1_imag = -GSL_IMAG(twiddle1[k - 1]);
      w2_real = GSL_REAL(twiddle2[k - 1]);
      w2_imag = -GSL_IMAG(twiddle2[k - 1]);
      w3_real = GSL_REAL(twiddle3[k - 1]);
      w3_imag = -GSL_IMAG(twiddle3[k - 1]);

      for (k1 = 0; k1 < q; k1++)
	{
	  const size_t from0 = k1 * product_1 + 2 * k - 1;
	  const size_t from1 = from0 + m;
	  const size_t from2 = from1 + m;
	  const size_t from3 = from2 + m;
	  
	  const double f0_real = VECTOR(in,istride,from0);
	  const double f0_imag = VECTOR(in,istride,from0 + 1);
	  const double f1_real = VECTOR(in,istride,from1);
	  const double f1_imag = VECTOR(in,istride,from1 + 1);
	  const double f2_real = VECTOR(in,istride,from2);
	  const double f2_imag = VECTOR(in,istride,from2 + 1);
	  const double f3_real = VECTOR(in,istride,from3);
	  const double f3_imag = VECTOR(in,istride,from3 + 1);
	  
	  const double z0_real = f0_real;
	  const double z0_imag = f0_imag;
	  const double z1_real = w1_real * f1_real - w1_imag * f1_imag;
	  const double z1_imag = w1_real * f1_imag + w1_imag * f1_real;
	  const double z2_real = w2_real * f2_real - w2_imag * f2_imag;
	  const double z2_imag = w2_real * f2_imag + w2_imag * f2_real;
	  const double z3_real = w3_real * f3_real - w3_imag * f3_imag;
	  const double z3_imag = w3_real * f3_imag + w3_imag * f3_real;

	  /* compute x = W(4) z */

	  /* t1 = z0 + z2 */
	  const double t1_real = z0_real + z2_real;
	  const double t1_imag = z0_imag + z2_imag;
	  
	  /* t2 = z1 + z3 */
	  const double t2_real = z1_real + z3_real;
	  const double t2_imag = z1_imag + z3_imag;
	  
	  /* t3 = z0 - z2 */
	  const double t3_real = z0_real - z2_real;
	  const double t3_imag = z0_imag - z2_imag;
	  
	  /* t4 = - (z1 - z3) */
	  const double t4_real = -(z1_real - z3_real);
	  const double t4_imag = -(z1_imag - z3_imag);
	  
	  /* x0 = t1 + t2 */
	  const double x0_real = t1_real + t2_real;
	  const double x0_imag = t1_imag + t2_imag;
	  
	  /* x1 = t3 + i t4 */
	  const double x1_real = t3_real - t4_imag;
	  const double x1_imag = t3_imag + t4_real;
	  
	  /* x2 = t1 - t2 */
	  const double x2_real = t1_real - t2_real;
	  const double x2_imag = t1_imag - t2_imag;
	  
	  /* x3 = t3 - i t4 */
	  const double x3_real = t3_real + t4_imag;
	  const double x3_imag = t3_imag - t4_real;

	  const size_t to0 = k1 * product + 2 * k - 1;
	  const size_t to1 = to0 + 2 * product_1;
	  const size_t to2 = 2 * product_1 - 2 * k + k1 * product - 1;
	  const size_t to3 = to2 + 2 * product_1;
	  
	  VECTOR(out,ostride,to0) = x0_real;
	  VECTOR(out,ostride,to0 + 1) = x0_imag;
	  
	  VECTOR(out,ostride,to1) = x1_real;
	  VECTOR(out,ostride,to1 + 1) = x1_imag;
	  
	  VECTOR(out,ostride,to3) = x2_real;
	  VECTOR(out,ostride,to3 + 1) = -x2_imag;
	  
	  VECTOR(out,ostride,to2) = x3_real;
	  VECTOR(out,ostride,to2 + 1) = -x3_imag;
	}
    }

  if (product_1 % 2 == 1)
    return;

  for (k1 = 0; k1 < q; k1++)
    {
      const size_t from0 = k1 * product_1 + product_1 - 1;
      const size_t from1 = from0 + m;
      const size_t from2 = from1 + m;
      const size_t from3 = from2 + m;
      
      const double x0 = VECTOR(in,istride,from0);
      const double x1 = VECTOR(in,istride,from1);
      const double x2 = VECTOR(in,istride,from2);
      const double x3 = VECTOR(in,istride,from3);
      
      const double t1 = (1.0 / sqrt (2.0)) * (x1 - x3);
      const double t2 = (1.0 / sqrt (2.0)) * (x1 + x3);
      
      const size_t to0 = k1 * product + 2 * k - 1;
      const size_t to1 = to0 + 2 * product_1;
      
      VECTOR(out,ostride,to0) = x0 + t1;
      VECTOR(out,ostride,to0 + 1) = -x2 - t2;
      
      VECTOR(out,ostride,to1) = x0 - t1;
      VECTOR(out,ostride,to1 + 1) = x2 - t2;
    }
  return;
}
