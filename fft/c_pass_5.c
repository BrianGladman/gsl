#include <config.h>
#include <math.h>

#include <gsl_complex.h>
#include <gsl_fft_complex.h>

#include "fft_complex.h"

int
gsl_fft_complex_pass_5 (const gsl_complex from[],
			gsl_complex to[],
			const gsl_fft_direction sign,
			const size_t product,
			const size_t n,
			const gsl_complex twiddle1[],
			const gsl_complex twiddle2[],
			const gsl_complex twiddle3[],
			const gsl_complex twiddle4[])
{
  size_t i = 0, j = 0;
  size_t k, k1;

  const size_t factor = 5;
  const size_t m = n / factor;
  const size_t q = n / product;
  const size_t product_1 = product / factor;
  const size_t jump = (factor - 1) * product_1;

  const double sin_2pi_by_5 = sin (2.0 * M_PI / 5.0);
  const double sin_2pi_by_10 = sin (2.0 * M_PI / 10.0);

  for (k = 0; k < q; k++)
    {

      double w1_real, w1_imag, w2_real, w2_imag, w3_real, w3_imag, w4_real,
        w4_imag;

      if (k == 0)
	{
	  w1_real = 1.0;
	  w1_imag = 0.0;
	  w2_real = 1.0;
	  w2_imag = 0.0;
	  w3_real = 1.0;
	  w3_imag = 0.0;
	  w4_real = 1.0;
	  w4_imag = 0.0;
	}
      else
	{
	  if (sign == forward)
	    {
	      /* forward tranform */
	      w1_real = twiddle1[k - 1].real;
	      w1_imag = twiddle1[k - 1].imag;
	      w2_real = twiddle2[k - 1].real;
	      w2_imag = twiddle2[k - 1].imag;
	      w3_real = twiddle3[k - 1].real;
	      w3_imag = twiddle3[k - 1].imag;
	      w4_real = twiddle4[k - 1].real;
	      w4_imag = twiddle4[k - 1].imag;
	    }
	  else
	    {
	      /* backward tranform: w -> conjugate(w) */
	      w1_real = twiddle1[k - 1].real;
	      w1_imag = -twiddle1[k - 1].imag;
	      w2_real = twiddle2[k - 1].real;
	      w2_imag = -twiddle2[k - 1].imag;
	      w3_real = twiddle3[k - 1].real;
	      w3_imag = -twiddle3[k - 1].imag;
	      w4_real = twiddle4[k - 1].real;
	      w4_imag = -twiddle4[k - 1].imag;
	    }
	}

      for (k1 = 0; k1 < product_1; k1++)
	{

	  gsl_complex z0, z1, z2, z3, z4;
	  double x0_real, x0_imag, x1_real, x1_imag, x2_real, x2_imag,
	    x3_real, x3_imag, x4_real, x4_imag;

	  {
	    const size_t from0 = i;
	    const size_t from1 = from0 + m;
	    const size_t from2 = from1 + m;
	    const size_t from3 = from2 + m;
	    const size_t from4 = from3 + m;

	    z0 = from[from0];
	    z1 = from[from1];
	    z2 = from[from2];
	    z3 = from[from3];
	    z4 = from[from4];
	  }

	  /* compute x = W(5) z */
	  {
	    /* t1 = z1 + z4 */
	    const double t1_real = z1.real + z4.real;
	    const double t1_imag = z1.imag + z4.imag;

	    /* t2 = z2 + z3 */
	    const double t2_real = z2.real + z3.real;
	    const double t2_imag = z2.imag + z3.imag;

	    /* t3 = z1 - z4 */
	    const double t3_real = z1.real - z4.real;
	    const double t3_imag = z1.imag - z4.imag;

	    /* t4 = z2 - z3 */
	    const double t4_real = z2.real - z3.real;
	    const double t4_imag = z2.imag - z3.imag;

	    /* t5 = t1 + t2 */
	    const double t5_real = t1_real + t2_real;
	    const double t5_imag = t1_imag + t2_imag;

	    /* t6 = (sqrt(5)/4)(t1 - t2) */
	    const double t6_real = (sqrt (5.0) / 4.0) * (t1_real - t2_real);
	    const double t6_imag = (sqrt (5.0) / 4.0) * (t1_imag - t2_imag);

	    /* t7 = z0 - ((t5)/4) */
	    const double t7_real = z0.real - t5_real / 4.0;
	    const double t7_imag = z0.imag - t5_imag / 4.0;

	    /* t8 = t7 + t6 */
	    const double t8_real = t7_real + t6_real;
	    const double t8_imag = t7_imag + t6_imag;

	    /* t9 = t7 - t6 */
	    const double t9_real = t7_real - t6_real;
	    const double t9_imag = t7_imag - t6_imag;

	    /* t10 = sin(2 pi/5) t3 + sin(2 pi/10) t4 */
	    const double t10_real = ((int) sign) * (sin_2pi_by_5 * t3_real +
						    sin_2pi_by_10 * t4_real);
	    const double t10_imag = ((int) sign) * (sin_2pi_by_5 * t3_imag +
						    sin_2pi_by_10 * t4_imag);

	    /* t11 = sin(2 pi/10) t3 - sin(2 pi/5) t4 */
	    const double t11_real = ((int) sign) * (sin_2pi_by_10 * t3_real -
						    sin_2pi_by_5 * t4_real);
	    const double t11_imag = ((int) sign) * (sin_2pi_by_10 * t3_imag -
						    sin_2pi_by_5 * t4_imag);

	    /* x0 = z0 + t5 */
	    x0_real = z0.real + t5_real;
	    x0_imag = z0.imag + t5_imag;

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

	  /* apply twiddle factors */
	  {
	    const size_t to0 = j;
	    const size_t to1 = to0 + product_1;
	    const size_t to2 = to1 + product_1;
	    const size_t to3 = to2 + product_1;
	    const size_t to4 = to3 + product_1;

	    /* to0 = 1 * x0 */
	    to[to0].real = x0_real;
	    to[to0].imag = x0_imag;

	    /* to1 = w1 * x1 */
	    to[to1].real = w1_real * x1_real - w1_imag * x1_imag;
	    to[to1].imag = w1_real * x1_imag + w1_imag * x1_real;

	    /* to2 = w2 * x2 */
	    to[to2].real = w2_real * x2_real - w2_imag * x2_imag;
	    to[to2].imag = w2_real * x2_imag + w2_imag * x2_real;

	    /* to3 = w3 * x3 */
	    to[to3].real = w3_real * x3_real - w3_imag * x3_imag;
	    to[to3].imag = w3_real * x3_imag + w3_imag * x3_real;

	    /* to4 = w4 * x4 */
	    to[to4].real = w4_real * x4_real - w4_imag * x4_imag;
	    to[to4].imag = w4_real * x4_imag + w4_imag * x4_real;
	  }

	  i++;
	  j++;
	}
      j += jump;
    }
  return 0;
}
