#include <math.h>
#include <gsl_complex.h>
#include <gsl_fft_complex.h>

int
gsl_fft_complex_pass_n (complex from[],
			complex to[],
			const gsl_fft_direction sign,
			const unsigned int factor,
			const unsigned int product,
			const unsigned int n,
			const complex twiddle[])
{
  unsigned int i = 0, j = 0;
  unsigned int k, k1;

  const unsigned int m = n / factor;
  const unsigned int q = n / product;
  const unsigned int product_1 = product / factor;
  const unsigned int jump = (factor - 1) * product_1;

  double w_real, w_imag;

  int e, e1, e2;

  for (k = 0; k < q; k++)
    {
      for (k1 = 0; k1 < product_1; k1++)
	{
	  to[i] = from[i];
	  for (e = 1; e < (factor - 1) / 2 + 1; e++)
	    {
	      const unsigned int j = i + e * m;
	      const unsigned int jc = i + (factor - e) * m;
	      to[j].real = from[j].real + from[jc].real;
	      to[j].imag = from[j].imag + from[jc].imag;
	      to[jc].real = from[j].real - from[jc].real;
	      to[jc].imag = from[j].imag - from[jc].imag;
	    }

	  for (e = 0; e < factor; e++)
	    {
	      double sum_real = to[i].real;
	      double sum_imag = to[i].imag;
	      for (e1 = 1; e1 < (factor - 1) / 2 + 1; e1++)
		{
		  complex xp = to[i + e1 * m];
		  complex xm = to[i + (factor - e1) * m];
		  double w_real = cos (2 * M_PI * (double) e * (double) e1 / (double) factor);
		  double w_imag = sin (2 * (int) sign * M_PI * (double) e * (double) e1 / (double) factor);
		  sum_real += w_real * xp.real - w_imag * xm.imag;
		  sum_imag += w_real * xp.imag + w_imag * xm.real;
		}
	      from[i + e * m].real = sum_real;
	      from[i + e * m].imag = sum_imag;
	    }

	  i++;
	}
    }

  i = 0;
  j = 0;

  for (k = 0; k < q; k++)
    {
      for (k1 = 0; k1 < product_1; k1++)
	{
	  for (e1 = 0; e1 < factor; e1++)
	    {
	      double w_real = cos (2 * M_PI * (double) k * (double) e1 / (double) (factor * q));
	      double w_imag = sin (2 * (int) sign * M_PI * (double) k * (double) e1 / (double) (factor * q));
	      double x_real = from[i + e1 * m].real;
	      double x_imag = from[i + e1 * m].imag;

	      to[j + e1 * product_1].real = w_real * x_real - w_imag * x_imag;
	      to[j + e1 * product_1].imag = w_real * x_imag + w_imag * x_real;
	    }
	  i++;
	  j++;
	}
      j += jump;
    }

  return 0;
}
