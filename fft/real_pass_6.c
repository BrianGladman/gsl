#include <math.h>
#include <gsl_complex.h>
#include <gsl_fft_real.h>

#include <fft_real.h>

int
gsl_fft_real_pass_6 (const double from[], double to[],
		     const unsigned int product,
		     const unsigned int n,
		     const complex twiddle1[],
		     const complex twiddle2[],
		     const complex twiddle3[],
		     const complex twiddle4[],
		     const complex twiddle5[])
{

  unsigned int i, j, k, k1, jump;
  unsigned int factor, q, m, product_1;
  unsigned int from0, from1, from2, from3, from4;
  unsigned int to0, to1, to2, to3, to4;

  double tau = sqrt (3.0) / 2.0;

  i = 0;
  j = 0;

  factor = 6;
  m = n / factor;
  q = n / product;
  product_1 = product / factor;
  jump = (factor - 1) * product_1;

  for (k1 = 0; k1 < q; k1++)
    {
      double t1_real, t2_real, t3_real, t4_real, t5_real, t6_real, t7_real,
        t8_real, t9_real, t10_real, t11_real;
      double z0_real, z1_real, z2_real, z3_real, z4_real;
      complex x0, x1, x2;

      from0 = k1 * product_1;
      from1 = from0 + m;
      from2 = from1 + m;
      from3 = from2 + m;
      from4 = from3 + m;
      from5 = from4 + m;

      z0_real = from[from0];
      z1_real = from[from1];
      z2_real = from[from2];
      z3_real = from[from3];
      z4_real = from[from4];
      z5_real = from[from5];

      /* compute x = W(6) z */
      /* W(6) is a combination of sums and differences of W(3) acting
         on the even and odd elements of z */

      /* ta1 = z2 + z4 */
      ta1.real = z2.real + z4.real;

      /* ta2 = z0 - ta1/2 */
      ta2.real = z0.real - ta1.real / 2;

      /* ta3 = (+/-) sin(pi/3)*(z2 - z4) */
      ta3.real = tau * (z2.real - z4.real);

      /* a0 = z0 + ta1 */
      a0.real = z0.real + ta1.real;

      /* a1 = ta2 + i ta3 */
      a1.real = ta2.real - ta3.imag;

      /* a2 = ta2 - i ta3 */
      a2.real = ta2.real + ta3.imag;

      /* tb1 = z5 + z1 */
      tb1.real = z5.real + z1.real;

      /* tb2 = z3 - tb1/2 */
      tb2.real = z3.real - tb1.real / 2;

      /* tb3 = -sin(pi/3)*(z5 - z1) */
      tb3.real = -tau * (z5.real - z1.real);

      /* b0 = z3 + tb1 */
      b0.real = z3.real + tb1.real;

      /* b1 = tb2 + i tb3 */
      b1.real = tb2.real;
      b1.imag = tb3.real;

      /* b2 = tb2 - i tb3 */
      b2.real = tb2.real;
      b2.imag = -tb3.real;

      /* x0 = a0 + b0 */
      x0.real = a0.real + b0.real;
      x0.imag = a0.imag + b0.imag;

      /* x4 = a1 + b1 */
      x4.real = a1.real + b1.real;
      x4.imag = a1.imag + b1.imag;

      /* x2 = a2 + b2 */
      x2.real = a2.real + b2.real;
      x2.imag = a2.imag + b2.imag;

      /* x3 = a0 - b0 */
      x3.real = a0.real - b0.real;
      x3.imag = a0.imag - b0.imag;

      /* x1 = a1 - b1 */
      x1.real = a1.real - b1.real;
      x1.imag = a1.imag - b1.imag;

      /* x5 = a2 - b2 */
      x5.real = a2.real - b2.real;
      x5.imag = a2.imag - b2.imag;

      to0 = product * k1;
      to1 = to0 + 2 * product_1 - 1;
      to2 = to1 + 2 * product_1;

      to[to0] = x0.real;
      to[to1] = x1.real;
      to[to1 + 1] = x1.imag;
      to[to2] = x2.real;
      to[to2 + 1] = x2.imag;
      to[to3] = x3.real;
      to[to3 + 1] = x3.imag;

    }

  if (product_1 == 1)
    return 0;

  for (k = 1; k < (product_1 + 1) / 2; k++)
    {
      complex w1, w2, w3, w4;
      w1.real = twiddle1[k - 1].real;
      w1.imag = -twiddle1[k - 1].imag;
      w2.real = twiddle2[k - 1].real;
      w2.imag = -twiddle2[k - 1].imag;
      w3.real = twiddle3[k - 1].real;
      w3.imag = -twiddle3[k - 1].imag;
      w4.real = twiddle4[k - 1].real;
      w4.imag = -twiddle4[k - 1].imag;

      for (k1 = 0; k1 < q; k1++)
	{
	  complex t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11;
	  complex z0, z1, z2, z3, z4;
	  complex x0, x1, x2, x3, x4;

	  from0 = k1 * product_1 + 2 * k - 1;
	  from1 = from0 + m;
	  from2 = from1 + m;
	  from3 = from2 + m;
	  from4 = from3 + m;

	  z0.real = from[from0];
	  z0.imag = from[from0 + 1];
	  z1.real = w1.real * from[from1] - w1.imag * from[from1 + 1];
	  z1.imag = w1.real * from[from1 + 1] + w1.imag * from[from1];
	  z2.real = w2.real * from[from2] - w2.imag * from[from2 + 1];
	  z2.imag = w2.real * from[from2 + 1] + w2.imag * from[from2];
	  z3.real = w3.real * from[from3] - w3.imag * from[from3 + 1];
	  z3.imag = w3.real * from[from3 + 1] + w3.imag * from[from3];
	  z4.real = w4.real * from[from4] - w4.imag * from[from4 + 1];
	  z4.imag = w4.real * from[from4 + 1] + w4.imag * from[from4];

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

	  /* FIXME signs */
	  /* t10 = sin(2 pi/5) t3 + sin(2 pi/10) t4 */
	  t10.real = -sina * t3.real - sinb * t4.real;
	  t10.imag = -sina * t3.imag - sinb * t4.imag;

	  /* FIXME signes */
	  /* t11 = sin(2 pi/10) t3 - sin(2 pi/5) t4 */
	  t11.real = -sinb * t3.real + sina * t4.real;
	  t11.imag = -sinb * t3.imag + sina * t4.imag;

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

	  to0 = k1 * product + 2 * k - 1;
	  to1 = to0 + 2 * product_1;
	  to2 = to1 + 2 * product_1;
	  to3 = 2 * product_1 - 2 * k + k1 * product - 1;
	  to4 = to3 + 2 * product_1;

	  /* to0 = 1 * x0 */
	  to[to0] = x0.real;
	  to[to0 + 1] = x0.imag;

	  /* to1 = 1 * x1 */
	  to[to1] = x1.real;
	  to[to1 + 1] = x1.imag;

	  /* to2 = 1 * x2 */
	  to[to2] = x2.real;
	  to[to2 + 1] = x2.imag;

	  to[to3] = x4.real;
	  to[to3 + 1] = -x4.imag;

	  to[to4] = x3.real;
	  to[to4 + 1] = -x3.imag;

	}
    }

  if (product_1 % 2 == 1)
    return 0;

  printf ("BARF q = %d, product_1 = %d\n", q, product_1);

  for (k1 = 0; k1 < q; k1++)
    {
      double x0, x1, x2, x3, x4;
      double t1, t2, t3, t4, t5, t6, t7;

      from0 = k1 * product_1 + product_1 - 1;
      from1 = from0 + m;
      from2 = from1 + m;
      from3 = from2 + m;
      from4 = from3 + m;

      x0 = from[from0];
      x1 = from[from1];
      x2 = from[from2];
      x3 = from[from3];
      x4 = from[from4];

      t1 = x1 - x4;
      t2 = x1 + x4;
      t3 = x2 - x3;
      t4 = x2 + x3;
      t5 = t1 - t3;
      t6 = x0 + t5 / 4.0;
      t7 = (sqrt (5.0) / 4.0) * (t1 + t3);

      to0 = k1 * product + product_1 - 1;
      to1 = to0 + 2 * product_1;
      to2 = to1 + 2 * product_1;
      to3 = 2 * product_1 - product_1 + k1 * product - 1;
      to4 = to3 + 2 * product_1;

      to[to0] = t6 + t7;
      to[to0 + 1] = -sinb * t2 - sina * t4;

      to[to1] = t6 - t7;
      to[to1 + 1] = -sina * t2 + sinb * t4;

      to[to2] = x0 - t5;
    }

  return 0;
}
