#include <config.h>
#include <math.h>

#include <gsl_complex.h>
#include <gsl_fft_real.h>

#include "fft_real.h"

int
gsl_fft_real_pass_4 (const double from[],
		     double to[],
		     const size_t product,
		     const size_t n,
		     const gsl_complex twiddle1[],
		     const gsl_complex twiddle2[],
		     const gsl_complex twiddle3[])
{
  size_t k, k1;

  const size_t factor = 4;
  const size_t m = n / factor;
  const size_t q = n / product;
  const size_t product_1 = product / factor;

  for (k1 = 0; k1 < q; k1++)
    {
      double z0_real, z1_real, z2_real, z3_real;
      double x0_real, x1_real, x1_imag, x2_real;

      {
	const size_t from0 = k1 * product_1;
	const size_t from1 = from0 + m;
	const size_t from2 = from1 + m;
	const size_t from3 = from2 + m;

	z0_real = from[from0];
	z1_real = from[from1];
	z2_real = from[from2];
	z3_real = from[from3];
      }

      /* compute x = W(4) z */
      {
	/* t1 = z0 + z2 */
	const double t1_real = z0_real + z2_real;

	/* t2 = z1 + z3 */
	const double t2_real = z1_real + z3_real;

	/* t3 = z0 - z2 */
	const double t3_real = z0_real - z2_real;

	/* t4 = - (z1 - z3) */
	const double t4_real = -(z1_real - z3_real);

	/* x0 = t1 + t2 */
	x0_real = t1_real + t2_real;

	/* x1 = t3 + i t4 */
	x1_real = t3_real;
	x1_imag = t4_real;

	/* x2 = t1 - t2 */
	x2_real = t1_real - t2_real;
      }

      {
	const size_t to0 = product * k1;
	const size_t to1 = to0 + 2 * product_1 - 1;
	const size_t to2 = to1 + 2 * product_1;

	to[to0] = x0_real;
	to[to1] = x1_real;
	to[to1 + 1] = x1_imag;
	to[to2] = x2_real;
      }
    }

  if (product_1 == 1)
    return 0;

  for (k = 1; k < (product_1 + 1) / 2; k++)
    {
      double w1_real, w1_imag, w2_real, w2_imag, w3_real, w3_imag;
      w1_real = twiddle1[k - 1].real;
      w1_imag = -twiddle1[k - 1].imag;
      w2_real = twiddle2[k - 1].real;
      w2_imag = -twiddle2[k - 1].imag;
      w3_real = twiddle3[k - 1].real;
      w3_imag = -twiddle3[k - 1].imag;

      for (k1 = 0; k1 < q; k1++)
	{
	  double z0_real, z0_imag, z1_real, z1_imag, z2_real, z2_imag,
	    z3_real, z3_imag;
	  double x0_real, x0_imag, x1_real, x1_imag, x2_real, x2_imag,
	    x3_real, x3_imag;

	  {
	    const size_t from0 = k1 * product_1 + 2 * k - 1;
	    const size_t from1 = from0 + m;
	    const size_t from2 = from1 + m;
	    const size_t from3 = from2 + m;

	    const double f0_real = from[from0];
	    const double f0_imag = from[from0 + 1];
	    const double f1_real = from[from1];
	    const double f1_imag = from[from1 + 1];
	    const double f2_real = from[from2];
	    const double f2_imag = from[from2 + 1];
	    const double f3_real = from[from3];
	    const double f3_imag = from[from3 + 1];

	    z0_real = f0_real;
	    z0_imag = f0_imag;
	    z1_real = w1_real * f1_real - w1_imag * f1_imag;
	    z1_imag = w1_real * f1_imag + w1_imag * f1_real;
	    z2_real = w2_real * f2_real - w2_imag * f2_imag;
	    z2_imag = w2_real * f2_imag + w2_imag * f2_real;
	    z3_real = w3_real * f3_real - w3_imag * f3_imag;
	    z3_imag = w3_real * f3_imag + w3_imag * f3_real;
	  }

	  /* compute x = W(4) z */
	  {
	    /* t1 = z0 + z2 */
	    const double t1_real = z0_real + z2_real;
	    const double t1_imag = z0_imag + z2_imag;

	    /* t2 = z1 + z3 */
	    const double t2_real = z1_real + z3_real;
	    const double t2_imag = z1_imag + z3_imag;

	    /* t3 = z0 - z2 */
	    const double t3_real = z0_real - z2_real;
	    const double t3_imag = z0_imag - z2_imag;

	    /* t4 = - (z1 - z3) */
	    const double t4_real = -(z1_real - z3_real);
	    const double t4_imag = -(z1_imag - z3_imag);

	    /* x0 = t1 + t2 */
	    x0_real = t1_real + t2_real;
	    x0_imag = t1_imag + t2_imag;

	    /* x1 = t3 + i t4 */
	    x1_real = t3_real - t4_imag;
	    x1_imag = t3_imag + t4_real;

	    /* x2 = t1 - t2 */
	    x2_real = t1_real - t2_real;
	    x2_imag = t1_imag - t2_imag;

	    /* x3 = t3 - i t4 */
	    x3_real = t3_real + t4_imag;
	    x3_imag = t3_imag - t4_real;
	  }

	  {
	    const size_t to0 = k1 * product + 2 * k - 1;
	    const size_t to1 = to0 + 2 * product_1;
	    const size_t to2 = 2 * product_1 - 2 * k + k1 * product - 1;
	    const size_t to3 = to2 + 2 * product_1;

	    to[to0] = x0_real;
	    to[to0 + 1] = x0_imag;

	    to[to1] = x1_real;
	    to[to1 + 1] = x1_imag;

	    to[to3] = x2_real;
	    to[to3 + 1] = -x2_imag;

	    to[to2] = x3_real;
	    to[to2 + 1] = -x3_imag;
	  }
	}
    }

  if (product_1 % 2 == 1)
    return 0;

  for (k1 = 0; k1 < q; k1++)
    {
      double t1, t2;
      double x0, x1, x2, x3;

      {
	const size_t from0 = k1 * product_1 + product_1 - 1;
	const size_t from1 = from0 + m;
	const size_t from2 = from1 + m;
	const size_t from3 = from2 + m;

	x0 = from[from0];
	x1 = from[from1];
	x2 = from[from2];
	x3 = from[from3];
      }

      t1 = (1.0 / sqrt (2.0)) * (x1 - x3);
      t2 = (1.0 / sqrt (2.0)) * (x1 + x3);

      {
	const size_t to0 = k1 * product + 2 * k - 1;
	const size_t to1 = to0 + 2 * product_1;

	to[to0] = x0 + t1;
	to[to0 + 1] = -x2 - t2;

	to[to1] = x0 - t1;
	to[to1 + 1] = x2 - t2;
      }
    }
  return 0;
}
