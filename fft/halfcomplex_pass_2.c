#include <math.h>
#include <gsl_complex.h>
#include <gsl_fft_halfcomplex.h>

#include <fft_halfcomplex.h>

int
gsl_fft_halfcomplex_pass_2 (const double from[],
			    double to[],
			    const unsigned int product,
			    const unsigned int n,
			    const gsl_complex twiddle[])
{

  unsigned int i, j, k, k1, jump;
  unsigned int factor, q, m, product_1;
  unsigned int from0, from1;
  unsigned int to0, to1;
  gsl_complex x0, x1;
  gsl_complex w;
  gsl_complex z0, z1;
  double r0, r1, s0, s1;
  i = 0;
  j = 0;

  factor = 2;
  m = n / factor;
  q = n / product;
  product_1 = product / factor;
  jump = (factor - 1) * q;

  for (k1 = 0; k1 < product_1; k1++)
    {
      from0 = 2 * k1 * q;
      from1 = from0 + 2 * q - 1;

      r0 = from[from0];
      r1 = from[from1];

      s0 = r0 + r1;
      s1 = r0 - r1;

      to0 = q * k1;
      to1 = to0 + m;

      to[to0] = s0;
      to[to1] = s1;

    }

  if (q == 1)
    return 0;

  for (k = 1; k < (q + 1) / 2; k++)
    {

      w.real = twiddle[k - 1].real;
      w.imag = twiddle[k - 1].imag;

      for (k1 = 0; k1 < product_1; k1++)
	{
	  from0 = 2 * k1 * q + 2 * k - 1;
	  from1 = 2 * k1 * q - 2 * k + 2 * q - 1;

	  z0.real = from[from0];
	  z0.imag = from[from0 + 1];

	  z1.real = from[from1];
	  z1.imag = from[from1 + 1];

	  /* compute x = W(2) z */

	  /* x0 = z0 + z1 */
	  x0.real = z0.real + z1.real;
	  x0.imag = z0.imag - z1.imag;

	  /* x1 = z0 - z1 */
	  x1.real = z0.real - z1.real;
	  x1.imag = z0.imag + z1.imag;

	  to0 = k1 * q + 2 * k - 1;
	  to1 = to0 + m;

	  to[to0] = x0.real;
	  to[to0 + 1] = x0.imag;

	  to[to1] = w.real * x1.real - w.imag * x1.imag;
	  to[to1 + 1] = w.imag * x1.real + w.real * x1.imag;

	}
    }

  if (q % 2 == 1)
    return 0;

  for (k1 = 0; k1 < product_1; k1++)
    {
      from0 = 2 * k1 * q + q - 1;
      to0 = k1 * q + q - 1;
      to1 = to0 + m;

      to[to0] = 2 * from[from0];
      to[to1] = -2 * from[from0 + 1];
    }
  return 0;
}
