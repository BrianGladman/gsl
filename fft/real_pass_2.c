#include <math.h>
#include <gsl_complex.h>
#include <gsl_fft_real.h>

#include <fft_real.h>

int
gsl_fft_real_pass_2 (const double from[],
		     double to[],
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
      double r0, r1, s0, s1;
      {
	const size_t from0 = k1 * product_1;
	const size_t from1 = from0 + m;

	r0 = from[from0];
	r1 = from[from1];
      }

      s0 = r0 + r1;
      s1 = r0 - r1;

      {
	const size_t to0 = product * k1;
	const size_t to1 = to0 + product - 1;

	to[to0] = s0;
	to[to1] = s1;
      }

    }

  if (product_1 == 1)
    return 0;

  for (k = 1; k < (product_1 + 1) / 2; k++)
    {

      /* forward real transform: w -> conjugate(w) */
      const double w_real = twiddle[k - 1].real;
      const double w_imag = -twiddle[k - 1].imag;

      for (k1 = 0; k1 < q; k1++)
	{
	  double f0_real, f0_imag, f1_real, f1_imag, z0_real, z0_imag,
	    z1_real, z1_imag, x0_real, x0_imag, x1_real, x1_imag;

	  {
	    const size_t from0 = k1 * product_1 + 2 * k - 1;
	    const size_t from1 = from0 + m;

	    f0_real = from[from0];
	    f0_imag = from[from0 + 1];

	    f1_real = from[from1];
	    f1_imag = from[from1 + 1];
	  }

	  z0_real = f0_real;
	  z0_imag = f0_imag;

	  z1_real = w_real * f1_real - w_imag * f1_imag;
	  z1_imag = w_real * f1_imag + w_imag * f1_real;

	  /* compute x = W(2) z */

	  /* x0 = z0 + z1 */
	  x0_real = z0_real + z1_real;
	  x0_imag = z0_imag + z1_imag;

	  /* x1 = z0 - z1 */
	  x1_real = z0_real - z1_real;
	  x1_imag = z0_imag - z1_imag;

	  {
	    const size_t to0 = k1 * product + 2 * k - 1;
	    const size_t to1 = k1 * product + product - 2 * k - 1;

	    to[to0] = x0_real;
	    to[to0 + 1] = x0_imag;

	    /* stored in conjugate location */
	    to[to1] = x1_real;
	    to[to1 + 1] = -x1_imag;
	  }
	}
    }

  if (product_1 % 2 == 1)
    return 0;

  for (k1 = 0; k1 < q; k1++)
    {
      const size_t from0 = k1 * product_1 + product_1 - 1;
      const size_t from1 = from0 + m;
      const size_t to0 = k1 * product + product_1 - 1;

      to[to0] = from[from0];
      to[to0 + 1] = -from[from1];
    }
  return 0;
}
