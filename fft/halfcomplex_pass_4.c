#include <math.h>
#include <gsl_complex.h>
#include <gsl_fft_halfcomplex.h>

#include <fft_halfcomplex.h>

int
gsl_fft_halfcomplex_pass_4 (const double from[],
			    double to[],
			    const size_t product,
			    const size_t n,
			    const gsl_complex twiddle1[],
			    const gsl_complex twiddle2[],
			    const gsl_complex twiddle3[])
{

  size_t i, j, k, k1, jump;
  size_t factor, q, m, product_1;
  size_t from0, from1, from2, from3;
  size_t to0, to1, to2, to3;

  i = 0;
  j = 0;

  factor = 4;
  m = n / factor;
  q = n / product;
  product_1 = product / factor;
  jump = (factor - 1) * q;

  for (k1 = 0; k1 < product_1; k1++)
    {
      double z0_real, z1_real, z1_imag, z2_real;
      double t1_real, t2_real, t3_real, t4_imag;

      from0 = 4 * k1 * q;
      from1 = from0 + 2 * q - 1;
      from2 = from1 + 2 * q;

      z0_real = from[from0];
      z1_real = from[from1];
      z1_imag = from[from1 + 1];
      z2_real = from[from2];

      t1_real = z0_real + z2_real;
      t2_real = 2 * z1_real;
      t3_real = z0_real - z2_real;
      t4_imag = 2 * z1_imag;

      to0 = q * k1;
      to1 = to0 + m;
      to2 = to1 + m;
      to3 = to2 + m;

      to[to0] = t1_real + t2_real;
      to[to1] = t3_real - t4_imag;
      to[to2] = t1_real - t2_real;
      to[to3] = t3_real + t4_imag;
    }

  if (q == 1)
    return 0;

  for (k = 1; k < (q + 1) / 2; k++)
    {
      gsl_complex w1, w2, w3;
      w1 = twiddle1[k - 1];
      w2 = twiddle2[k - 1];
      w3 = twiddle3[k - 1];

      for (k1 = 0; k1 < product_1; k1++)
	{
	  gsl_complex t1, t2, t3, t4;
	  gsl_complex z0, z1, z2, z3;
	  gsl_complex x0, x1, x2, x3;

	  from0 = 4 * k1 * q + 2 * k - 1;
	  from1 = from0 + 2 * q;
	  from2 = 4 * k1 * q - 2 * k + 2 * q - 1;
	  from3 = from2 + 2 * q;

	  z0.real = from[from0];
	  z0.imag = from[from0 + 1];

	  z1.real = from[from1];
	  z1.imag = from[from1 + 1];

	  z2.real = from[from3];
	  z2.imag = -from[from3 + 1];

	  z3.real = from[from2];
	  z3.imag = -from[from2 + 1];

	  /* compute x = W(4) z */

	  /* t1 = z0 + z2 */
	  t1.real = z0.real + z2.real;
	  t1.imag = z0.imag + z2.imag;

	  /* t2 = z1 + z3 */
	  t2.real = z1.real + z3.real;
	  t2.imag = z1.imag + z3.imag;

	  /* t3 = z0 - z2 */
	  t3.real = z0.real - z2.real;
	  t3.imag = z0.imag - z2.imag;

	  /* t4 = (z1 - z3) */
	  t4.real = (z1.real - z3.real);
	  t4.imag = (z1.imag - z3.imag);

	  /* x0 = t1 + t2 */
	  x0.real = t1.real + t2.real;
	  x0.imag = t1.imag + t2.imag;

	  /* x1 = t3 + i t4 */
	  x1.real = t3.real - t4.imag;
	  x1.imag = t3.imag + t4.real;

	  /* x2 = t1 - t2 */
	  x2.real = t1.real - t2.real;
	  x2.imag = t1.imag - t2.imag;


	  /* x3 = t3 - i t4 */
	  x3.real = t3.real + t4.imag;
	  x3.imag = t3.imag - t4.real;

	  to0 = k1 * q + 2 * k - 1;
	  to1 = to0 + m;
	  to2 = to1 + m;
	  to3 = to2 + m;

	  to[to0] = x0.real;
	  to[to0 + 1] = x0.imag;

	  to[to1] = w1.real * x1.real - w1.imag * x1.imag;
	  to[to1 + 1] = w1.imag * x1.real + w1.real * x1.imag;

	  to[to2] = w2.real * x2.real - w2.imag * x2.imag;
	  to[to2 + 1] = w2.imag * x2.real + w2.real * x2.imag;

	  /* to3 = w3 * x3 */
	  to[to3] = w3.real * x3.real - w3.imag * x3.imag;
	  to[to3 + 1] = w3.real * x3.imag + w3.imag * x3.real;

	}
    }

  if (q % 2 == 1)
    return 0;

  for (k1 = 0; k1 < product_1; k1++)
    {
      double z0_real, z0_imag, z1_real, z1_imag;
      double x0_real, x1_real, x2_real, x3_real;
      double t1_real, t2_real;

      from0 = 4 * k1 * q + q - 1;
      from1 = from0 + 2 * q;

      z0_real = from[from0];
      z0_imag = from[from0 + 1];

      z1_real = from[from1];
      z1_imag = from[from1 + 1];

      t1_real = sqrt (2.0) * (z0_imag + z1_imag);
      t2_real = sqrt (2.0) * (z0_real - z1_real);

      x0_real = 2 * (z0_real + z1_real);
      x1_real = t2_real - t1_real;
      x2_real = 2 * (z1_imag - z0_imag);
      x3_real = -(t2_real + t1_real);

      to0 = k1 * q + q - 1;
      to1 = to0 + m;
      to2 = to1 + m;
      to3 = to2 + m;

      to[to0] = x0_real;
      to[to1] = x1_real;
      to[to2] = x2_real;
      to[to3] = x3_real;
    }
  return 0;
}
