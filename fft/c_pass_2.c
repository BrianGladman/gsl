#include <config.h>
#include <math.h>

#include <gsl_complex.h>
#include <gsl_fft_complex.h>

#include "fft_complex.h"

int
fft_complex_pass_2 (const double in[],
		    const size_t istride,
		    double out[],
		    const size_t ostride,
		    const gsl_fft_direction sign,
		    const size_t product,
		    const size_t n,
		    const gsl_complex twiddle[])
{
  size_t i = 0, j = 0;
  size_t k, k1;

  const size_t factor = 2;
  const size_t m = n / factor;
  const size_t q = n / product;
  const size_t product_1 = product / factor;
  const size_t jump = (factor - 1) * product_1;

  for (k = 0; k < q; k++)
    {
      double w_real, w_imag;

      if (k == 0)
	{
	  w_real = 1.0;
	  w_imag = 0.0;
	}
      else
	{
	  if (sign == forward)
	    {
	      /* forward tranform */
	      w_real = GSL_REAL(twiddle[k - 1]);
	      w_imag = GSL_IMAG(twiddle[k - 1]);
	    }
	  else
	    {
	      /* backward tranform: w -> conjugate(w) */
	      w_real = GSL_REAL(twiddle[k - 1]);
	      w_imag = -GSL_IMAG(twiddle[k - 1]);
	    }
	}

      for (k1 = 0; k1 < product_1; k1++)
	{
	  const double z0_real = REAL(in,istride,i);
	  const double z0_imag = IMAG(in,istride,i);

	  const double z1_real = REAL(in,istride,i+m);
	  const double z1_imag = IMAG(in,istride,i+m);

	  /* compute x = W(2) z */

	  /* x0 = z0 + z1 */
	  const double x0_real = z0_real + z1_real;
	  const double x0_imag = z0_imag + z1_imag;

	  /* x1 = z0 - z1 */
	  const double x1_real = z0_real - z1_real;
	  const double x1_imag = z0_imag - z1_imag;

	  /* apply twiddle factors */
	  
	  /* out0 = 1 * x0 */
	  REAL(out,ostride,j) = x0_real;
	  IMAG(out,ostride,j) = x0_imag;
	  
	  /* out1 = w * x1 */
	  REAL(out,ostride,j+product_1) = w_real * x1_real - w_imag * x1_imag;
	  IMAG(out,ostride,j+product_1) = w_real * x1_imag + w_imag * x1_real;
	  
	  i++;
	  j++;
	}
      j += jump;
    }
  return 0;
}
