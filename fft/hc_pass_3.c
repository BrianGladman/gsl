#include <config.h>
#include <math.h>

#include <gsl_complex.h>
#include <gsl_fft_halfcomplex.h>

#include "fft_halfcomplex.h"

int
gsl_fft_halfcomplex_pass_3 (const double from[],
			    double to[],
			    const size_t product,
			    const size_t n,
			    const gsl_complex twiddle1[],
			    const gsl_complex twiddle2[])
{
  size_t i, j, k, k1, jump;
  size_t factor, q, m, product_1;

  double tau = sqrt (3.0) / 2.0;

  i = 0;
  j = 0;

  factor = 3;
  m = n / factor;
  q = n / product;
  product_1 = product / factor;
  jump = (factor - 1) * q;

  for (k1 = 0; k1 < product_1; k1++)
    {
      const size_t from0 = 3 * k1 * q;
      const size_t from1 = from0 + 2 * q - 1;

      const double z0_real = from[from0];
      const double z1_real = from[from1];
      const double z1_imag = from[from1 + 1];

      const double t1_real = 2 * z1_real;
      const double t2_real = z0_real - z1_real;
      const double t3_imag = 2 * tau * z1_imag;

      const size_t to0 = q * k1;
      const size_t to1 = to0 + m;
      const size_t to2 = to1 + m;

      to[to0] = z0_real + t1_real;
      to[to1] = t2_real - t3_imag;
      to[to2] = t2_real + t3_imag;

    }

  if (q == 1)
    return 0;

  for (k = 1; k < (q + 1) / 2; k++)
    {
      const double w1_real = GSL_REAL(twiddle1[k - 1]);
      const double w1_imag = GSL_IMAG(twiddle1[k - 1]);
      const double w2_real = GSL_REAL(twiddle2[k - 1]);
      const double w2_imag = GSL_IMAG(twiddle2[k - 1]);

      for (k1 = 0; k1 < product_1; k1++)
	{
	  const size_t from0 = 3 * k1 * q + 2 * k - 1;
	  const size_t from1 = from0 + 2 * q;
	  const size_t from2 = 3 * k1 * q - 2 * k + 2 * q - 1;

	  const double z0_real = from[from0];
	  const double z0_imag = from[from0 + 1];

	  const double z1_real = from[from1];
	  const double z1_imag = from[from1 + 1];

	  const double z2_real = from[from2];
	  const double z2_imag = -from[from2 + 1];

	  /* compute x = W(3) z */

	  /* t1 = z1 + z2 */
	  const double t1_real = z1_real + z2_real;
	  const double t1_imag = z1_imag + z2_imag;

	  /* t2 = z0 - t1/2 */
	  const double t2_real = z0_real - t1_real / 2.0;
	  const double t2_imag = z0_imag - t1_imag / 2.0;

	  /* t3 = sin(pi/3)*(z1 - z2) */
	  const double t3_real = tau * (z1_real - z2_real);
	  const double t3_imag = tau * (z1_imag - z2_imag);

	  /* x0 = z0 + t1 */
	  const double x0_real = z0_real + t1_real;
	  const double x0_imag = z0_imag + t1_imag;

	  /* x1 = t2 + i t3 */
	  const double x1_real = t2_real - t3_imag;
	  const double x1_imag = t2_imag + t3_real;

	  /* x2 = t2 - i t3 */
	  const double x2_real = t2_real + t3_imag;
	  const double x2_imag = t2_imag - t3_real;

	  const size_t to0 = k1 * q + 2 * k - 1;
	  const size_t to1 = to0 + m;
	  const size_t to2 = to1 + m;

	  to[to0] = x0_real;
	  to[to0 + 1] = x0_imag;

	  to[to1] = w1_real * x1_real - w1_imag * x1_imag;
	  to[to1 + 1] = w1_imag * x1_real + w1_real * x1_imag;

	  to[to2] = w2_real * x2_real - w2_imag * x2_imag;
	  to[to2 + 1] = w2_imag * x2_real + w2_real * x2_imag;

	}
    }

  if (q % 2 == 1)
    return 0;

  for (k1 = 0; k1 < product_1; k1++)
    {
      const size_t from0 = 3 * k1 * q + q - 1;
      const size_t from1 = from0 + 2 * q;

      const double z0_real = from[from0];
      const double z0_imag = from[from0 + 1];
      const double z1_real = from[from1];

      const double t1_real = z0_real - z1_real;
      const double t2_real = 2 * tau * z0_imag;

      const double x0_real = 2 * z0_real + z1_real;
      const double x1_real = t1_real - t2_real;
      const double x2_real = -t1_real - t2_real;

      const size_t to0 = k1 * q + q - 1;
      const size_t to1 = to0 + m;
      const size_t to2 = to1 + m;

      to[to0] = x0_real;
      to[to1] = x1_real;
      to[to2] = x2_real;
    }
  return 0;
}
