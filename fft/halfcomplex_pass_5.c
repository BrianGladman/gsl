#include <math.h>
#include <gsl_complex.h>
#include <gsl_fft_halfcomplex.h>

#include <fft_halfcomplex.h>

int
gsl_fft_halfcomplex_pass_5 (const double from[],
			    double to[],
			    const unsigned int product,
			    const unsigned int n,
			    const gsl_complex twiddle1[],
			    const gsl_complex twiddle2[],
			    const gsl_complex twiddle3[],
			    const gsl_complex twiddle4[])
{

  unsigned int i, j, k, k1, jump;
  unsigned int factor, q, m, product_1;
  unsigned int from0, from1, from2, from3, from4;
  unsigned int to0, to1, to2, to3, to4;

  const double sina = sin (2.0 * M_PI / 5.0);
  const double sinb = sin (2.0 * M_PI / 10.0);

  i = 0;
  j = 0;

  factor = 5;
  m = n / factor;
  q = n / product;
  product_1 = product / factor;
  jump = (factor - 1) * q;

  for (k1 = 0; k1 < product_1; k1++)
    {
      double z0_real, z1_real, z1_imag, z2_real, z2_imag;
      double t1_real, t2_real, t3_real, t4_real, t5_real, t6_imag, t7_imag;
      double x0_real, x1_real, x2_real, x3_real, x4_real;

      from0 = 5 * k1 * q;
      from1 = from0 + 2 * q - 1;
      from2 = from1 + 2 * q;

      z0_real = from[from0];
      z1_real = from[from1];
      z1_imag = from[from1 + 1];
      z2_real = from[from2];
      z2_imag = from[from2 + 1];


      t1_real = 2 * (z1_real + z2_real);
      t2_real = 2 * (sqrt (5.0) / 4.0) * (z1_real - z2_real);
      t3_real = z0_real - t1_real / 4.0;
      t4_real = t2_real + t3_real;
      t5_real = -t2_real + t3_real;
      t6_imag = 2 * (sina * z1_imag + sinb * z2_imag);
      t7_imag = 2 * (sinb * z1_imag - sina * z2_imag);

      x0_real = z0_real + t1_real;
      x1_real = t4_real - t6_imag;
      x2_real = t5_real - t7_imag;
      x3_real = t5_real + t7_imag;
      x4_real = t4_real + t6_imag;

      to0 = q * k1;
      to1 = to0 + m;
      to2 = to1 + m;
      to3 = to2 + m;
      to4 = to3 + m;

      to[to0] = x0_real;
      to[to1] = x1_real;
      to[to2] = x2_real;
      to[to3] = x3_real;
      to[to4] = x4_real;


    }

  if (q == 1)
    return 0;

  for (k = 1; k < (q + 1) / 2; k++)
    {
      gsl_complex w1, w2, w3, w4;
      w1 = twiddle1[k - 1];
      w2 = twiddle2[k - 1];
      w3 = twiddle3[k - 1];
      w4 = twiddle4[k - 1];

      for (k1 = 0; k1 < product_1; k1++)
	{
	  gsl_complex t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11;
	  gsl_complex z0, z1, z2, z3, z4;
	  gsl_complex x0, x1, x2, x3, x4;

	  from0 = 5 * k1 * q + 2 * k - 1;
	  from1 = from0 + 2 * q;
	  from2 = from1 + 2 * q;
	  from3 = 5 * k1 * q - 2 * k + 2 * q - 1;
	  from4 = from3 + 2 * q;

	  z0.real = from[from0];
	  z0.imag = from[from0 + 1];

	  z1.real = from[from1];
	  z1.imag = from[from1 + 1];

	  z2.real = from[from2];
	  z2.imag = from[from2 + 1];

	  z3.real = from[from4];
	  z3.imag = -from[from4 + 1];

	  z4.real = from[from3];
	  z4.imag = -from[from3 + 1];

	  /* compute x = W(5) z */

	  /* t1 = z1 + z4 */
	  t1.real = z1.real + z4.real;
	  t1.imag = z1.imag + z4.imag;

	  /* t2 = z2 + z3 */
	  t2.real = z2.real + z3.real;
	  t2.imag = z2.imag + z3.imag;

	  /* t3 = z1 - z4 */
	  t3.real = z1.real - z4.real;
	  t3.imag = z1.imag - z4.imag;

	  /* t4 = z2 - z3 */
	  t4.real = z2.real - z3.real;
	  t4.imag = z2.imag - z3.imag;

	  /* t5 = t1 + t2 */
	  t5.real = t1.real + t2.real;
	  t5.imag = t1.imag + t2.imag;

	  /* t6 = (sqrt(5)/4)(t1 - t2) */
	  t6.real = (sqrt (5.0) / 4.0) * (t1.real - t2.real);
	  t6.imag = (sqrt (5.0) / 4.0) * (t1.imag - t2.imag);

	  /* t7 = z0 - ((t5)/4) */
	  t7.real = z0.real - t5.real / 4.0;
	  t7.imag = z0.imag - t5.imag / 4.0;

	  /* t8 = t7 + t6 */
	  t8.real = t7.real + t6.real;
	  t8.imag = t7.imag + t6.imag;

	  /* t9 = t7 - t6 */
	  t9.real = t7.real - t6.real;
	  t9.imag = t7.imag - t6.imag;

	  /* t10 = sin(2 pi/5) t3 + sin(2 pi/10) t4 */
	  t10.real = sina * t3.real + sinb * t4.real;
	  t10.imag = sina * t3.imag + sinb * t4.imag;

	  /* t11 = sin(2 pi/10) t3 - sin(2 pi/5) t4 */
	  t11.real = sinb * t3.real - sina * t4.real;
	  t11.imag = sinb * t3.imag - sina * t4.imag;

	  /* x0 = z0 + t5 */
	  x0.real = z0.real + t5.real;
	  x0.imag = z0.imag + t5.imag;

	  /* x1 = t8 + i t10 */
	  x1.real = t8.real - t10.imag;
	  x1.imag = t8.imag + t10.real;

	  /* x2 = t9 + i t11 */
	  x2.real = t9.real - t11.imag;
	  x2.imag = t9.imag + t11.real;

	  /* x3 = t9 - i t11 */
	  x3.real = t9.real + t11.imag;
	  x3.imag = t9.imag - t11.real;

	  /* x4 = t8 - i t10 */
	  x4.real = t8.real + t10.imag;
	  x4.imag = t8.imag - t10.real;

	  to0 = k1 * q + 2 * k - 1;
	  to1 = to0 + m;
	  to2 = to1 + m;
	  to3 = to2 + m;
	  to4 = to3 + m;

	  /* apply twiddle factors */

	  /* to0 = 1 * x0 */
	  to[to0] = x0.real;
	  to[to0 + 1] = x0.imag;

	  /* to1 = w1 * x1 */
	  to[to1] = w1.real * x1.real - w1.imag * x1.imag;
	  to[to1 + 1] = w1.real * x1.imag + w1.imag * x1.real;

	  /* to2 = w2 * x2 */
	  to[to2] = w2.real * x2.real - w2.imag * x2.imag;
	  to[to2 + 1] = w2.real * x2.imag + w2.imag * x2.real;

	  /* to3 = w3 * x3 */
	  to[to3] = w3.real * x3.real - w3.imag * x3.imag;
	  to[to3 + 1] = w3.real * x3.imag + w3.imag * x3.real;

	  /* to4 = w4 * x4 */
	  to[to4] = w4.real * x4.real - w4.imag * x4.imag;
	  to[to4 + 1] = w4.real * x4.imag + w4.imag * x4.real;

	}
    }

  if (q % 2 == 1)
    return 0;

  for (k1 = 0; k1 < product_1; k1++)
    {
      double z0_real, z0_imag, z1_real, z1_imag, z2_real;
      double x0_real, x1_real, x2_real, x3_real, x4_real;
      double t1_real, t2_real, t3_real, t4_real, t5_real, t6_real, t7_real;

      from0 = 5 * k1 * q + q - 1;
      from1 = from0 + 2 * q;
      from2 = from1 + 2 * q;

      z0_real = 2 * from[from0];
      z0_imag = 2 * from[from0 + 1];

      z1_real = 2 * from[from1];
      z1_imag = 2 * from[from1 + 1];

      z2_real = from[from2];

      t1_real = z0_real + z1_real;
      t2_real = (t1_real / 4.0) - z2_real;
      t3_real = (sqrt (5.0) / 4.0) * (z0_real - z1_real);
      t4_real = sinb * z0_imag + sina * z1_imag;
      t5_real = sina * z0_imag - sinb * z1_imag;
      t6_real = t3_real + t2_real;
      t7_real = t3_real - t2_real;

      x0_real = t1_real + z2_real;
      x1_real = t6_real - t4_real;
      x2_real = t7_real - t5_real;
      x3_real = -t7_real - t5_real;
      x4_real = -t6_real - t4_real;

      to0 = k1 * q + q - 1;
      to1 = to0 + m;
      to2 = to1 + m;
      to3 = to2 + m;
      to4 = to3 + m;

      to[to0] = x0_real;
      to[to1] = x1_real;
      to[to2] = x2_real;
      to[to3] = x3_real;
      to[to4] = x4_real;
    }
  return 0;
}
