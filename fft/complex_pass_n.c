#include <math.h>
#include <gsl_complex.h>
#include <gsl_fft_complex.h>

int
gsl_fft_complex_pass_n (const complex from[],
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

  for (k = 0; k < q; k++)
    {
      for (k1 = 0; k1 < product_1; k1++)
	{
	  /* compute x = W(factor) z */

	  unsigned int e1;

	  for (e1 = 0; e1 < factor; e1++)
	    {
	      unsigned int e2;
	      unsigned int idx = 0;
	      const unsigned int idx_step = e1 * q;
	      double sum_real = 0.0, sum_imag = 0.0;

	      for (e2 = 0; e2 < factor; e2++)
		{
		  double t_real, t_imag;
		  const unsigned int from0 = e2 * m + i;
		  double f_real = from[from0].real;
		  double f_imag = from[from0].imag;

		  if (idx == 0)
		    {
		      t_real = 1.0;
		      t_imag = 0.0;
		    }
		  else
		    {
		      /* make use of twiddle vector to compute
		         elements of W matrix */
		      if (sign == forward)
			{
			  t_real = twiddle[idx - 1].real;
			  t_imag = twiddle[idx - 1].imag;
			}
		      else
			{
			  t_real = twiddle[idx - 1].real;
			  t_imag = -twiddle[idx - 1].imag;
			}
		    }

		  sum_real += t_real * f_real - t_imag * f_imag;
		  sum_imag += t_real * f_imag + t_imag * f_real;
		  idx += idx_step;
		  idx %= factor * q;
		}

	      /* apply twiddle factor */

	      if (e1 == 0 || k == 0)
		{
		  w_real = 1.0;
		  w_imag = 0.0;
		}
	      else
		{
		  const unsigned int w_idx = (e1 - 1) * q + k - 1;
		  if (sign == forward)
		    {
		      w_real = twiddle[w_idx].real;
		      w_imag = twiddle[w_idx].imag;
		    }
		  else
		    {
		      w_real = twiddle[w_idx].real;
		      w_imag = -twiddle[w_idx].imag;
		    }
		};

	      {
		const unsigned int to_idx = product_1 * e1 + j;

		to[to_idx].real = w_real * sum_real - w_imag * sum_imag;
		to[to_idx].imag = w_real * sum_imag + w_imag * sum_real;
	      }

	    }
	  i++;
	  j++;
	}
      j += jump;
    }
  return 0;
}
