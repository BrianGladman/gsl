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

	 
	  /* e = 0 */
	  {
	    double sum_real = to[i].real ;
	    double sum_imag = to[i].imag ;
	    for (e1 = 1; e1 < (factor - 1) / 2 + 1; e1++)
	      {
		sum_real += to[i + e1*m].real ;
		sum_imag += to[i + e1*m].imag ;
	      }
	    from[i].real = sum_real ;
	    from[i].imag = sum_imag ;
	  }

	  for (e = 1; e < (factor-1)/2 + 1; e++)
	    {
	      double sum_real = 0 ;
	      double sum_imag = 0 ;
	      double sumc_real = 0 ;
	      double sumc_imag = 0 ;
	      unsigned int idx = e*q ;
	      const unsigned int idx_step = e * q ;
	      for (e1 = 1; e1 < (factor - 1) / 2 + 1; e1++)
		{
		  complex xp = to[i + e1 * m];
		  complex xm = to[i + (factor - e1) * m];
		  double w_real, w_imag ;

		  if (idx == 0) {
		    w_real = 1 ;
		    w_imag = 0 ;
		  } else {
		    if (sign == forward) {
		      w_real = twiddle[idx - 1].real ;
		      w_imag = twiddle[idx - 1].imag ;
		    } else {
		      w_real = twiddle[idx - 1].real ;
		      w_imag = -twiddle[idx - 1].imag ;
		    }
		  }

		  sum_real += w_real * xp.real - w_imag * xm.imag;
		  sum_imag += w_real * xp.imag + w_imag * xm.real;
		  sumc_real += w_real * xp.real + w_imag * xm.imag;
		  sumc_imag += w_real * xp.imag - w_imag * xm.real;
		  idx += idx_step ;
		  idx %= factor * q ;
		}
	      from[i + e * m].real = to[i].real + sum_real;
	      from[i + e * m].imag = to[i].imag + sum_imag;
	      from[i + (factor-e) * m].real = to[i].real + sumc_real;
	      from[i + (factor-e) * m].imag = to[i].imag + sumc_imag;
	    }

	  i++;
	}
    }

  i = 0;
  j = 0;

  /* k = 0 */
  for (k1 = 0; k1 < product_1; k1++)
    {
      
      to[j].real = from[i].real;
      to[j].imag = from[i].imag;
      
      for (e1 = 1; e1 < factor; e1++)
	{
	  double x_real = from[i + e1 * m].real;
	  double x_imag = from[i + e1 * m].imag;
	  to[j + e1 * product_1].real = x_real;
	  to[j + e1 * product_1].imag = x_imag;
	}
      i++;
      j++;
    }
  j += jump;

  for (k = 1; k < q; k++)
    {
      for (k1 = 0; k1 < product_1; k1++)
	{

	  to[j].real = from[i].real;
	  to[j].imag = from[i].imag;

	  for (e1 = 1; e1 < factor; e1++)
	    {
	      double x_real = from[i + e1 * m].real;
	      double x_imag = from[i + e1 * m].imag;

	      double wtr, wti ;
	      if (k == 0) {
		w_real = 1 ;
		w_imag = 0 ;
	      } else {
		if (sign == forward) {
		  w_real = twiddle[(e1-1)*q + k-1].real ;
		  w_imag = twiddle[(e1-1)*q + k-1].imag ;
		} else {
		  w_real = twiddle[(e1-1)*q + k-1].real ;
		  w_imag = -twiddle[(e1-1)*q + k-1].imag ; 
		}
	      }

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

