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

  unsigned int e, e1;

  for (i = 0; i < m; i++)
    {
      to[i] = from[i];
    }

  for (e = 1; e < (factor - 1) / 2 + 1; e++)
    {
      for (i = 0; i < m; i++)
	{
	  const unsigned int idx = i + e * m;
	  const unsigned int idxc = i + (factor - e) * m;
	  to[idx].real = from[idx].real + from[idxc].real;
	  to[idx].imag = from[idx].imag + from[idxc].imag;
	  to[idxc].real = from[idx].real - from[idxc].real;
	  to[idxc].imag = from[idx].imag - from[idxc].imag;
	}
    }

  /* e = 0 */

  for (i=0 ; i<m; i++) 
    {
      from[i] = to[i] ;
    }

  for (e1 = 1; e1 < (factor - 1) / 2 + 1; e1++)
    {
      for (i = 0; i < m; i++)
	{
	  from[i].real += to[i + e1*m].real ;
	  from[i].imag += to[i + e1*m].imag ;
	}
    }

  for (e = 1; e < (factor-1)/2 + 1; e++)
    {
      unsigned int idx = e*q ;
      const unsigned int idx_step = e * q ;
      double w_real, w_imag ;

      const unsigned int em = e * m ;
      const unsigned int ecm = (factor - e) * m ;

      for (i = 0; i < m; i++) 
	{
	  from[i + em] = to[i];
	  from[i + ecm] = to[i];
	}

      for (e1 = 1; e1 < (factor - 1) / 2 + 1; e1++)
	{
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

	  for (i = 0; i < m; i++) 
	    {
	      complex xp = to[i + e1 * m];
	      complex xm = to[i + (factor - e1) *m];
	
	      const double ap = w_real * xp.real ;
	      const double am = w_imag * xm.imag ; 

	      double sum_real = ap - am;
	      double sumc_real = ap + am;

	      const double bp = w_real * xp.imag ;
	      const double bm = w_imag * xm.real ;

	      double sum_imag = bp + bm;
	      double sumc_imag = bp - bm;

	      from[i + em].real += sum_real;
	      from[i + em].imag += sum_imag;
	      from[i + ecm].real +=  sumc_real;
	      from[i + ecm].imag += sumc_imag;
	    }
	  idx += idx_step ;
	  idx %= factor * q ;
	}
    }

  i = 0;
  j = 0;

  /* k = 0 */
  for (k1 = 0; k1 < product_1; k1++)
    {
      to[k1] = from[k1];
    }

  for (e1 = 1; e1 < factor; e1++)
    {
      for (k1 = 0; k1 < product_1; k1++)
	{
	  to[k1 + e1 * product_1] = from[k1 + e1 * m] ;
	}
    }

  i = product_1 ;
  j = product ;

  for (k = 1; k < q; k++)
    {
      for (k1 = 0; k1 < product_1; k1++)
	{

	  to[j].real = from[i].real;
	  to[j].imag = from[i].imag;

	  i++;
	  j++;
	}
      j += jump;
    }

  i = product_1 ;
  j = product ;

  for (k = 1; k < q; k++)
    {
      for (k1 = 0; k1 < product_1; k1++)
	{
	  for (e1 = 1; e1 < factor; e1++)
	    {
	      double x_real = from[i + e1 * m].real;
	      double x_imag = from[i + e1 * m].imag;

	      double w_real, w_imag ;
	      if (sign == forward) {
		w_real = twiddle[(e1-1)*q + k-1].real ;
		w_imag = twiddle[(e1-1)*q + k-1].imag ;
	      } else {
		w_real = twiddle[(e1-1)*q + k-1].real ;
		w_imag = -twiddle[(e1-1)*q + k-1].imag ; 
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

