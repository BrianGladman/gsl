#include <config.h>
#include <math.h>

#include <gsl_complex.h>
#include <gsl_fft_complex.h>

#include "fft_complex.h"

int
gsl_fft_complex_pass_6 (const gsl_complex from[],
			gsl_complex to[],
			const gsl_fft_direction sign,
			const size_t product,
			const size_t n,
			const gsl_complex twiddle1[],
			const gsl_complex twiddle2[],
			const gsl_complex twiddle3[],
			const gsl_complex twiddle4[],
			const gsl_complex twiddle5[])
{

  size_t i = 0, j = 0;
  size_t k, k1;

  const size_t factor = 6;
  const size_t m = n / factor;
  const size_t q = n / product;
  const size_t product_1 = product / factor;
  const size_t jump = (factor - 1) * product_1;

  const double tau = sqrt (3.0) / 2.0;

  for (k = 0; k < q; k++)
    {
      double w1_real, w1_imag, w2_real, w2_imag, w3_real, w3_imag, w4_real,
        w4_imag, w5_real, w5_imag;

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
	  w5_real = 1.0;
	  w5_imag = 0.0;
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
	      w5_real = twiddle5[k - 1].real;
	      w5_imag = twiddle5[k - 1].imag;
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
	      w5_real = twiddle5[k - 1].real;
	      w5_imag = -twiddle5[k - 1].imag;
	    }
	}

      for (k1 = 0; k1 < product_1; k1++)
	{
	  gsl_complex z0, z1, z2, z3, z4, z5;
	  double x0_real, x0_imag, x1_real, x1_imag, x2_real, x2_imag,
	    x3_real, x3_imag, x4_real, x4_imag, x5_real, x5_imag;

	  {
	    const size_t from0 = i;
	    const size_t from1 = from0 + m;
	    const size_t from2 = from1 + m;
	    const size_t from3 = from2 + m;
	    const size_t from4 = from3 + m;
	    const size_t from5 = from4 + m;

	    z0 = from[from0];
	    z1 = from[from1];
	    z2 = from[from2];
	    z3 = from[from3];
	    z4 = from[from4];
	    z5 = from[from5];
	  }

	  /* compute x = W(6) z */

	  /* W(6) is a combination of sums and differences of W(3) acting
	     on the even and odd elements of z */
	  {
	    /* ta1 = z2 + z4 */
	    const double ta1_real = z2.real + z4.real;
	    const double ta1_imag = z2.imag + z4.imag;

	    /* ta2 = z0 - ta1/2 */
	    const double ta2_real = z0.real - ta1_real / 2;
	    const double ta2_imag = z0.imag - ta1_imag / 2;

	    /* ta3 = (+/-) sin(pi/3)*(z2 - z4) */
	    const double ta3_real = ((int) sign) * tau * (z2.real - z4.real);
	    const double ta3_imag = ((int) sign) * tau * (z2.imag - z4.imag);

	    /* a0 = z0 + ta1 */
	    const double a0_real = z0.real + ta1_real;
	    const double a0_imag = z0.imag + ta1_imag;

	    /* a1 = ta2 + i ta3 */
	    const double a1_real = ta2_real - ta3_imag;
	    const double a1_imag = ta2_imag + ta3_real;

	    /* a2 = ta2 - i ta3 */
	    const double a2_real = ta2_real + ta3_imag;
	    const double a2_imag = ta2_imag - ta3_real;

	    /* tb1 = z5 + z1 */
	    const double tb1_real = z5.real + z1.real;
	    const double tb1_imag = z5.imag + z1.imag;

	    /* tb2 = z3 - tb1/2 */
	    const double tb2_real = z3.real - tb1_real / 2;
	    const double tb2_imag = z3.imag - tb1_imag / 2;

	    /* tb3 = (+/-) sin(pi/3)*(z5 - z1) */
	    const double tb3_real = ((int) sign) * tau * (z5.real - z1.real);
	    const double tb3_imag = ((int) sign) * tau * (z5.imag - z1.imag);

	    /* b0 = z3 + tb1 */
	    const double b0_real = z3.real + tb1_real;
	    const double b0_imag = z3.imag + tb1_imag;

	    /* b1 = tb2 + i tb3 */
	    const double b1_real = tb2_real - tb3_imag;
	    const double b1_imag = tb2_imag + tb3_real;

	    /* b2 = tb2 - i tb3 */
	    const double b2_real = tb2_real + tb3_imag;
	    const double b2_imag = tb2_imag - tb3_real;

	    /* x0 = a0 + b0 */
	    x0_real = a0_real + b0_real;
	    x0_imag = a0_imag + b0_imag;

	    /* x4 = a1 + b1 */
	    x4_real = a1_real + b1_real;
	    x4_imag = a1_imag + b1_imag;

	    /* x2 = a2 + b2 */
	    x2_real = a2_real + b2_real;
	    x2_imag = a2_imag + b2_imag;

	    /* x3 = a0 - b0 */
	    x3_real = a0_real - b0_real;
	    x3_imag = a0_imag - b0_imag;

	    /* x1 = a1 - b1 */
	    x1_real = a1_real - b1_real;
	    x1_imag = a1_imag - b1_imag;

	    /* x5 = a2 - b2 */
	    x5_real = a2_real - b2_real;
	    x5_imag = a2_imag - b2_imag;
	  }

	  /* apply twiddle factors */
	  {
	    const size_t to0 = j;
	    const size_t to1 = to0 + product_1;
	    const size_t to2 = to1 + product_1;
	    const size_t to3 = to2 + product_1;
	    const size_t to4 = to3 + product_1;
	    const size_t to5 = to4 + product_1;

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

	    /* to5 = w5 * x5 */
	    to[to5].real = w5_real * x5_real - w5_imag * x5_imag;
	    to[to5].imag = w5_real * x5_imag + w5_imag * x5_real;
	  }
	  i++;
	  j++;
	}
      j += jump;
    }
  return 0;
}
