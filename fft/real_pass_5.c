#include <math.h>
#include <gsl_complex.h>
#include <gsl_fft_real.h>

int
gsl_fft_real_pass_5 (const double from[],
		     double to[],
		     const unsigned int product,
		     const unsigned int n,
		     const complex twiddle1[],
		     const complex twiddle2[],
		     const complex twiddle3[],
		     const complex twiddle4[])
{
  unsigned int k, k1;

  const unsigned int factor = 5;
  const unsigned int m = n / factor;
  const unsigned int q = n / product;
  const unsigned int product_1 = product / factor;

  const double sina = sin (2.0 * M_PI / 5.0);
  const double sinb = sin (2.0 * M_PI / 10.0);

  for (k1 = 0; k1 < q; k1++)
    {
      double t1_real, t2_real, t3_real, t4_real, t5_real, t6_real, t7_real,
        t8_real, t9_real, t10_real, t11_real;
      double z0_real, z1_real, z2_real, z3_real, z4_real;
      double x0_real, x1_real, x1_imag, x2_real, x2_imag;

      {
	const unsigned int from0 = k1 * product_1;
	const unsigned int from1 = from0 + m;
	const unsigned int from2 = from1 + m;
	const unsigned int from3 = from2 + m;
	const unsigned int from4 = from3 + m;

	z0_real = from[from0];
	z1_real = from[from1];
	z2_real = from[from2];
	z3_real = from[from3];
	z4_real = from[from4];
      }

      /* t1 = z1 + z4 */
      t1_real = z1_real + z4_real;

      /* t2 = z2 + z3 */
      t2_real = z2_real + z3_real;

      /* t3 = z1 - z4 */
      t3_real = z1_real - z4_real;

      /* t4 = z2 - z3 */
      t4_real = z2_real - z3_real;

      /* t5 = t1 + t2 */
      t5_real = t1_real + t2_real;

      /* t6 = (sqrt(5)/4)(t1 - t2) */
      t6_real = (sqrt (5.0) / 4.0) * (t1_real - t2_real);

      /* t7 = z0 - ((t5)/4) */
      t7_real = z0_real - t5_real / 4.0;

      /* t8 = t7 + t6 */
      t8_real = t7_real + t6_real;

      /* t9 = t7 - t6 */
      t9_real = t7_real - t6_real;

      /* t10 = -(sin(2 pi/5) t3 + sin(2 pi/10) t4 ) */
      t10_real = -sina * t3_real - sinb * t4_real;

      /* t11 = -(sin(2 pi/10) t3 - sin(2 pi/5) t4) */
      t11_real = -sinb * t3_real + sina * t4_real;

      /* x0 = z0 + t5 */
      x0_real = z0_real + t5_real;

      /* x1 = t8 + i t10 */
      x1_real = t8_real;
      x1_imag = t10_real;

      /* x2 = t9 + i t11 */
      x2_real = t9_real;
      x2_imag = t11_real;

      {
	const unsigned int to0 = product * k1;
	const unsigned int to1 = to0 + 2 * product_1 - 1;
	const unsigned int to2 = to1 + 2 * product_1;

	to[to0] = x0_real;
	to[to1] = x1_real;
	to[to1 + 1] = x1_imag;
	to[to2] = x2_real;
	to[to2 + 1] = x2_imag;
      }

    }

  if (product_1 == 1)
    return 0;

  for (k = 1; k < (product_1 + 1) / 2; k++)
    {
      const double w1_real = twiddle1[k - 1].real;
      const double w1_imag = -twiddle1[k - 1].imag;
      const double w2_real = twiddle2[k - 1].real;
      const double w2_imag = -twiddle2[k - 1].imag;
      const double w3_real = twiddle3[k - 1].real;
      const double w3_imag = -twiddle3[k - 1].imag;
      const double w4_real = twiddle4[k - 1].real;
      const double w4_imag = -twiddle4[k - 1].imag;

      for (k1 = 0; k1 < q; k1++)
	{
	  double x0_real, x0_imag, x1_real, x1_imag, x2_real, x2_imag,
	    x3_real, x3_imag, x4_real, x4_imag;
	  double z0_real, z0_imag, z1_real, z1_imag, z2_real, z2_imag,
	    z3_real, z3_imag, z4_real, z4_imag;
	  {
	    const unsigned int from0 = k1 * product_1 + 2 * k - 1;
	    const unsigned int from1 = from0 + m;
	    const unsigned int from2 = from1 + m;
	    const unsigned int from3 = from2 + m;
	    const unsigned int from4 = from3 + m;

	    const double f0_real = from[from0];
	    const double f0_imag = from[from0 + 1];
	    const double f1_real = from[from1];
	    const double f1_imag = from[from1 + 1];
	    const double f2_real = from[from2];
	    const double f2_imag = from[from2 + 1];
	    const double f3_real = from[from3];
	    const double f3_imag = from[from3 + 1];
	    const double f4_real = from[from4];
	    const double f4_imag = from[from4 + 1];

	    z0_real = f0_real;
	    z0_imag = f0_imag;
	    z1_real = w1_real * f1_real - w1_imag * f1_imag;
	    z1_imag = w1_real * f1_imag + w1_imag * f1_real;
	    z2_real = w2_real * f2_real - w2_imag * f2_imag;
	    z2_imag = w2_real * f2_imag + w2_imag * f2_real;
	    z3_real = w3_real * f3_real - w3_imag * f3_imag;
	    z3_imag = w3_real * f3_imag + w3_imag * f3_real;
	    z4_real = w4_real * f4_real - w4_imag * f4_imag;
	    z4_imag = w4_real * f4_imag + w4_imag * f4_real;
	  }

	  /* compute x = W(5) z */
	  {
	    /* t1 = z1 + z4 */
	    const double t1_real = z1_real + z4_real;
	    const double t1_imag = z1_imag + z4_imag;

	    /* t2 = z2 + z3 */
	    const double t2_real = z2_real + z3_real;
	    const double t2_imag = z2_imag + z3_imag;

	    /* t3 = z1 - z4 */
	    const double t3_real = z1_real - z4_real;
	    const double t3_imag = z1_imag - z4_imag;

	    /* t4 = z2 - z3 */
	    const double t4_real = z2_real - z3_real;
	    const double t4_imag = z2_imag - z3_imag;

	    /* t5 = t1 + t2 */
	    const double t5_real = t1_real + t2_real;
	    const double t5_imag = t1_imag + t2_imag;

	    /* t6 = (sqrt(5)/4)(t1 - t2) */
	    const double t6_real = (sqrt (5.0) / 4.0) * (t1_real - t2_real);
	    const double t6_imag = (sqrt (5.0) / 4.0) * (t1_imag - t2_imag);

	    /* t7 = z0 - ((t5)/4) */
	    const double t7_real = z0_real - t5_real / 4.0;
	    const double t7_imag = z0_imag - t5_imag / 4.0;

	    /* t8 = t7 + t6 */
	    const double t8_real = t7_real + t6_real;
	    const double t8_imag = t7_imag + t6_imag;

	    /* t9 = t7 - t6 */
	    const double t9_real = t7_real - t6_real;
	    const double t9_imag = t7_imag - t6_imag;

	    /* t10 = - (sin(2 pi/5) t3 + sin(2 pi/10) t4) */
	    const double t10_real = -sina * t3_real - sinb * t4_real;
	    const double t10_imag = -sina * t3_imag - sinb * t4_imag;

	    /* t11 = -(sin(2 pi/10) t3 - sin(2 pi/5) t4) */
	    const double t11_real = -sinb * t3_real + sina * t4_real;
	    const double t11_imag = -sinb * t3_imag + sina * t4_imag;

	    /* x0 = z0 + t5 */
	    x0_real = z0_real + t5_real;
	    x0_imag = z0_imag + t5_imag;

	    /* x1 = t8 + i t10 */
	    x1_real = t8_real - t10_imag;
	    x1_imag = t8_imag + t10_real;

	    /* x2 = t9 + i t11 */
	    x2_real = t9_real - t11_imag;
	    x2_imag = t9_imag + t11_real;

	    /* x3 = t9 - i t11 */
	    x3_real = t9_real + t11_imag;
	    x3_imag = t9_imag - t11_real;

	    /* x4 = t8 - i t10 */
	    x4_real = t8_real + t10_imag;
	    x4_imag = t8_imag - t10_real;
	  }

	  {
	    const unsigned int to0 = k1 * product + 2 * k - 1;
	    const unsigned int to1 = to0 + 2 * product_1;
	    const unsigned int to2 = to1 + 2 * product_1;
	    const unsigned int to3 = 2 * product_1 - 2 * k + k1 * product - 1;
	    const unsigned int to4 = to3 + 2 * product_1;

	    to[to0] = x0_real;
	    to[to0 + 1] = x0_imag;

	    to[to1] = x1_real;
	    to[to1 + 1] = x1_imag;

	    to[to2] = x2_real;
	    to[to2 + 1] = x2_imag;

	    to[to3] = x4_real;
	    to[to3 + 1] = -x4_imag;

	    to[to4] = x3_real;
	    to[to4 + 1] = -x3_imag;
	  }
	}
    }

  if (product_1 % 2 == 1)
    return 0;

  for (k1 = 0; k1 < q; k1++)
    {
      double z0_real, z1_real, z2_real, z3_real, z4_real;
      double t1, t2, t3, t4, t5, t6, t7;

      {
	const unsigned int from0 = k1 * product_1 + product_1 - 1;
	const unsigned int from1 = from0 + m;
	const unsigned int from2 = from1 + m;
	const unsigned int from3 = from2 + m;
	const unsigned int from4 = from3 + m;

	z0_real = from[from0];
	z1_real = from[from1];
	z2_real = from[from2];
	z3_real = from[from3];
	z4_real = from[from4];
      }

      t1 = z1_real - z4_real;
      t2 = z1_real + z4_real;
      t3 = z2_real - z3_real;
      t4 = z2_real + z3_real;
      t5 = t1 - t3;
      t6 = z0_real + t5 / 4.0;
      t7 = (sqrt (5.0) / 4.0) * (t1 + t3);

      {
	const unsigned int to0 = k1 * product + product_1 - 1;
	const unsigned int to1 = to0 + 2 * product_1;
	const unsigned int to2 = to1 + 2 * product_1;
	/* const unsigned int to3 = 2 * product_1 - product_1 + k1 * product - 1;
	   const unsigned int to4 = to3 + 2 * product_1; */

	to[to0] = t6 + t7;
	to[to0 + 1] = -sinb * t2 - sina * t4;

	to[to1] = t6 - t7;
	to[to1 + 1] = -sina * t2 + sinb * t4;

	to[to2] = z0_real - t5;
      }
    }

  return 0;
}
