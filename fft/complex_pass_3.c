#include <math.h>
#include <gsl_complex.h>
#include <gsl_fft_complex.h>

#include <fft_complex.h>

int
gsl_fft_complex_pass_3 (const gsl_complex from[],
			gsl_complex to[],
			const gsl_fft_direction sign,
			const unsigned int product,
			const unsigned int n,
			const gsl_complex * twiddle1,
			const gsl_complex * twiddle2)
{
  unsigned int i = 0, j = 0;
  unsigned int k, k1;

  const unsigned int factor = 3;
  const unsigned int m = n / factor;
  const unsigned int q = n / product;
  const unsigned int product_1 = product / factor;
  const unsigned int jump = (factor - 1) * product_1;

  const double tau = sqrt (3.0) / 2.0;

  for (k = 0; k < q; k++)
    {
      double w1_real, w1_imag, w2_real, w2_imag;

      if (k == 0)
	{
	  w1_real = 1.0;
	  w1_imag = 0.0;
	  w2_real = 1.0;
	  w2_imag = 0.0;
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
	    }
	  else
	    {
	      /* backward tranform: w -> conjugate(w) */
	      w1_real = twiddle1[k - 1].real;
	      w1_imag = -twiddle1[k - 1].imag;
	      w2_real = twiddle2[k - 1].real;
	      w2_imag = -twiddle2[k - 1].imag;
	    }
	}

      for (k1 = 0; k1 < product_1; k1++)
	{

	  gsl_complex z0, z1, z2;
	  double x0_real, x0_imag, x1_real, x1_imag, x2_real, x2_imag;

	  {
	    const unsigned int from0 = i;
	    const unsigned int from1 = from0 + m;
	    const unsigned int from2 = from1 + m;

	    z0 = from[from0];
	    z1 = from[from1];
	    z2 = from[from2];
	  }

	  /* compute x = W(3) z */
	  {
	    /* t1 = z1 + z2 */
	    const double t1_real = z1.real + z2.real;
	    const double t1_imag = z1.imag + z2.imag;

	    /* t2 = z0 - t1/2 */
	    const double t2_real = z0.real - t1_real / 2.0;
	    const double t2_imag = z0.imag - t1_imag / 2.0;

	    /* t3 = (+/-) sin(pi/3)*(z1 - z2) */
	    const double t3_real = ((int) sign) * tau * (z1.real - z2.real);
	    const double t3_imag = ((int) sign) * tau * (z1.imag - z2.imag);

	    /* x0 = z0 + t1 */
	    x0_real = z0.real + t1_real;
	    x0_imag = z0.imag + t1_imag;

	    /* x1 = t2 + i t3 */
	    x1_real = t2_real - t3_imag;
	    x1_imag = t2_imag + t3_real;

	    /* x2 = t2 - i t3 */
	    x2_real = t2_real + t3_imag;
	    x2_imag = t2_imag - t3_real;
	  }

	  /* apply twiddle factors */
	  {
	    const unsigned int to0 = j;
	    const unsigned int to1 = to0 + product_1;
	    const unsigned int to2 = to1 + product_1;

	    /* to0 = 1 * x0 */
	    to[to0].real = x0_real;
	    to[to0].imag = x0_imag;

	    /* to1 = w1 * x1 */
	    to[to1].real = w1_real * x1_real - w1_imag * x1_imag;
	    to[to1].imag = w1_real * x1_imag + w1_imag * x1_real;

	    /* to2 = w2 * x2 */
	    to[to2].real = w2_real * x2_real - w2_imag * x2_imag;
	    to[to2].imag = w2_real * x2_imag + w2_imag * x2_real;
	  }

	  i++;
	  j++;
	}
      j += jump;
    }
  return 0;
}
