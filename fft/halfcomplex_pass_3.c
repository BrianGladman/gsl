#include <math.h>
#include <gsl_complex.h>
#include <gsl_fft_halfcomplex.h>

#include <fft_halfcomplex.h>

int
gsl_fft_halfcomplex_pass_3 (const double from[],
			    double to[],
			    const unsigned int product,
			    const unsigned int n,
			    const gsl_complex twiddle1[],
			    const gsl_complex twiddle2[])
{

  unsigned int i, j, k, k1, jump;
  unsigned int factor, q, m, product_1;
  unsigned int from0, from1, from2;
  unsigned int to0, to1, to2;

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
      double z0_real, z1_real, z1_imag;
      double t1_real, t2_real, t3_imag;

      from0 = 3 * k1 * q;
      from1 = from0 + 2 * q - 1;

      z0_real = from[from0];
      z1_real = from[from1];
      z1_imag = from[from1 + 1];

      t1_real = 2 * z1_real;
      t2_real = z0_real - z1_real;
      t3_imag = 2 * tau * z1_imag;

      to0 = q * k1;
      to1 = to0 + m;
      to2 = to1 + m;

      to[to0] = z0_real + t1_real;
      to[to1] = t2_real - t3_imag;
      to[to2] = t2_real + t3_imag;

    }

  if (q == 1)
    return 0;

  for (k = 1; k < (q + 1) / 2; k++)
    {
      gsl_complex w1, w2;
      w1 = twiddle1[k - 1];
      w2 = twiddle2[k - 1];

      for (k1 = 0; k1 < product_1; k1++)
	{
	  gsl_complex t1, t2, t3;
	  gsl_complex z0, z1, z2;
	  gsl_complex x0, x1, x2;

	  from0 = 3 * k1 * q + 2 * k - 1;
	  from1 = from0 + 2 * q;
	  from2 = 3 * k1 * q - 2 * k + 2 * q - 1;

	  z0.real = from[from0];
	  z0.imag = from[from0 + 1];

	  z1.real = from[from1];
	  z1.imag = from[from1 + 1];

	  z2.real = from[from2];
	  z2.imag = -from[from2 + 1];

	  /* compute x = W(3) z */

	  /* t1 = z1 + z2 */
	  t1.real = z1.real + z2.real;
	  t1.imag = z1.imag + z2.imag;

	  /* t2 = z0 - t1/2 */
	  t2.real = z0.real - t1.real / 2.0;
	  t2.imag = z0.imag - t1.imag / 2.0;

	  /* t3 = sin(pi/3)*(z1 - z2) */
	  t3.real = tau * (z1.real - z2.real);
	  t3.imag = tau * (z1.imag - z2.imag);

	  /* x0 = z0 + t1 */
	  x0.real = z0.real + t1.real;
	  x0.imag = z0.imag + t1.imag;

	  /* x1 = t2 + i t3 */
	  x1.real = t2.real - t3.imag;
	  x1.imag = t2.imag + t3.real;

	  /* x2 = t2 - i t3 */
	  x2.real = t2.real + t3.imag;
	  x2.imag = t2.imag - t3.real;

	  to0 = k1 * q + 2 * k - 1;
	  to1 = to0 + m;
	  to2 = to1 + m;

	  to[to0] = x0.real;
	  to[to0 + 1] = x0.imag;

	  to[to1] = w1.real * x1.real - w1.imag * x1.imag;
	  to[to1 + 1] = w1.imag * x1.real + w1.real * x1.imag;

	  to[to2] = w2.real * x2.real - w2.imag * x2.imag;
	  to[to2 + 1] = w2.imag * x2.real + w2.real * x2.imag;

	}
    }

  if (q % 2 == 1)
    return 0;

  for (k1 = 0; k1 < product_1; k1++)
    {
      double z0_real, z0_imag, z1_real;
      double x0_real, x1_real, x2_real;
      double t1_real, t2_real;

      from0 = 3 * k1 * q + q - 1;
      from1 = from0 + 2 * q;

      z0_real = from[from0];
      z0_imag = from[from0 + 1];
      z1_real = from[from1];

      t1_real = z0_real - z1_real;
      t2_real = 2 * tau * z0_imag;

      x0_real = 2 * z0_real + z1_real;
      x1_real = t1_real - t2_real;
      x2_real = -t1_real - t2_real;

      to0 = k1 * q + q - 1;
      to1 = to0 + m;
      to2 = to1 + m;

      to[to0] = x0_real;
      to[to1] = x1_real;
      to[to2] = x2_real;
    }
  return 0;
}
