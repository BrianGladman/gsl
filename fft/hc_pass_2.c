#include <config.h>
#include <math.h>

#include <gsl_complex.h>
#include <gsl_fft_halfcomplex.h>

#include "fft_halfcomplex.h"

int
gsl_fft_halfcomplex_pass_2 (const double from[],
			    double to[],
			    const size_t product,
			    const size_t n,
			    const gsl_complex twiddle[])
{
  size_t i, j, k, k1, jump;
  size_t factor, q, m, product_1;
  i = 0;
  j = 0;

  factor = 2;
  m = n / factor;
  q = n / product;
  product_1 = product / factor;
  jump = (factor - 1) * q;

  for (k1 = 0; k1 < product_1; k1++)
    {
      const double r0 = from[2 * k1 * q];
      const double r1 = from[2 * k1 * q + 2 * q - 1];

      const double s0 = r0 + r1;
      const double s1 = r0 - r1;

      to[q * k1] = s0;
      to[q * k1 + m] = s1;
    }

  if (q == 1)
    return 0;

  for (k = 1; k < (q + 1) / 2; k++)
    {

      const double w_real = GSL_REAL(twiddle[k - 1]);
      const double w_imag = GSL_IMAG(twiddle[k - 1]);

      for (k1 = 0; k1 < product_1; k1++)
	{
	  const size_t from0 = 2 * k1 * q + 2 * k - 1;
	  const size_t from1 = 2 * k1 * q - 2 * k + 2 * q - 1;

	  const double z0_real = from[from0];
	  const double z0_imag = from[from0 + 1];

	  const double z1_real = from[from1];
	  const double z1_imag = from[from1 + 1];

	  /* compute x = W(2) z */

	  /* x0 = z0 + z1 */
	  const double x0_real = z0_real + z1_real;
	  const double x0_imag = z0_imag - z1_imag;

	  /* x1 = z0 - z1 */
	  const double x1_real = z0_real - z1_real;
	  const double x1_imag = z0_imag + z1_imag;

	  const size_t to0 = k1 * q + 2 * k - 1;
	  const size_t to1 = to0 + m;

	  to[to0] = x0_real;
	  to[to0 + 1] = x0_imag;

	  to[to1] = w_real * x1_real - w_imag * x1_imag;
	  to[to1 + 1] = w_imag * x1_real + w_real * x1_imag;

	}
    }

  if (q % 2 == 1)
    return 0;

  for (k1 = 0; k1 < product_1; k1++)
    {
      const size_t from0 = 2 * k1 * q + q - 1;
      const size_t to0 = k1 * q + q - 1;
      const size_t to1 = to0 + m;

      to[to0] = 2 * from[from0];
      to[to1] = -2 * from[from0 + 1];
    }
  return 0;
}
