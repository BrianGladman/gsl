#include <config.h>
#include <math.h>

#include <gsl_complex.h>
#include <gsl_fft_halfcomplex.h>

#include "fft_halfcomplex.h"

int
gsl_fft_halfcomplex_pass_n (const double in[],
			    const size_t istride,
			    double out[],
			    const size_t ostride,
			    const size_t factor,
			    const size_t product,
			    const size_t n,
			    const gsl_complex twiddle[])
{

  size_t k, k1;

  const size_t m = n / factor;
  const size_t q = n / product;
  const size_t product_1 = product / factor;

  size_t e1, e2;

  const double d_theta = 2.0 * M_PI / ((double) factor);
  const double cos_d_theta = cos (d_theta);
  const double sin_d_theta = sin (d_theta);

  for (k1 = 0; k1 < product_1; k1++)
    {
      /* compute z = W(factor) x, for x halfcomplex */

      double dw_real = 1.0, dw_imag = 0.0;

      for (e1 = 0; e1 < factor; e1++)
	{
	  double sum_real = 0.0;
	  double w_real = 1.0, w_imag = 0.0;

	  if (e1 > 0)
	    {
	      double tmp_real = dw_real * cos_d_theta - dw_imag * sin_d_theta;
	      double tmp_imag = dw_real * sin_d_theta + dw_imag * cos_d_theta;
	      dw_real = tmp_real;
	      dw_imag = tmp_imag;
	    }

	  for (e2 = 0; e2 <= factor - e2; e2++)
	    {
	      double z_real, z_imag;

	      if (e2 > 0)
		{
		  double tmp_real = dw_real * w_real - dw_imag * w_imag;
		  double tmp_imag = dw_real * w_imag + dw_imag * w_real;
		  w_real = tmp_real;
		  w_imag = tmp_imag;
		}

	      if (e2 == 0)
		{
		  size_t from_idx = factor * k1 * q;
		  z_real = VECTOR(in,istride,from_idx);
		  z_imag = 0.0;
		  sum_real += w_real * z_real - w_imag * z_imag;
		}
	      else if (e2 == factor - e2)
		{
		  size_t from_idx = factor * q * k1 + 2 * e2 * q - 1;
		  z_real = VECTOR(in,istride,from_idx);
		  z_imag = 0.0;
		  sum_real += w_real * z_real;
		}
	      else
		{
		  size_t from_idx = factor * q * k1 + 2 * e2 * q - 1;
		  z_real = VECTOR(in,istride,from_idx);
		  z_imag = VECTOR(in,istride,from_idx + 1);
		  sum_real += 2 * (w_real * z_real - w_imag * z_imag);
		}

	    }

	  {
	    const size_t to_idx = q * k1 + e1 * m;
	    VECTOR(out,ostride,to_idx) = sum_real;
	  }
	}
    }

  if (q == 1)
    return 0;

  for (k = 1; k < (q + 1) / 2; k++)
    {
      for (k1 = 0; k1 < product_1; k1++)
	{

	  double dw_real = 1.0, dw_imag = 0.0;

	  for (e1 = 0; e1 < factor; e1++)
	    {
	      double z_real, z_imag;
	      double sum_real = 0.0;
	      double sum_imag = 0.0;
	      double w_real = 1.0, w_imag = 0.0;

	      if (e1 > 0)
		{
		  double t_real = dw_real * cos_d_theta - dw_imag * sin_d_theta;
		  double t_imag = dw_real * sin_d_theta + dw_imag * cos_d_theta;
		  dw_real = t_real;
		  dw_imag = t_imag;
		}

	      for (e2 = 0; e2 < factor; e2++)
		{

		  if (e2 > 0)
		    {
		      double tmp_real = dw_real * w_real - dw_imag * w_imag;
		      double tmp_imag = dw_real * w_imag + dw_imag * w_real;
		      w_real = tmp_real;
		      w_imag = tmp_imag;
		    }

		  if (e2 < factor - e2)
		    {
		      const size_t from0 = factor * k1 * q + 2 * k + 2 * e2 * q - 1;
		      z_real = VECTOR(in,istride,from0);
		      z_imag = VECTOR(in,istride,from0 + 1);
		    }
		  else
		    {
		      const size_t from0 = factor * k1 * q - 2 * k + 2 * (factor - e2) * q - 1;
		      z_real = VECTOR(in,istride,from0);
		      z_imag = -VECTOR(in,istride,from0 + 1);
		    }

		  sum_real += w_real * z_real - w_imag * z_imag;
		  sum_imag += w_real * z_imag + w_imag * z_real;
		}

	      if (k == 0 || e1 == 0)
		{
		  w_real = 1.0;
		  w_imag = 0.0;
		}
	      else
		{
		  size_t tskip = (q + 1) / 2 - 1;
		  w_real = GSL_REAL(twiddle[k - 1 + tskip * (e1 - 1)]);
		  w_imag = GSL_IMAG(twiddle[k - 1 + tskip * (e1 - 1)]);
		}

	      {
		const size_t to0 = k1 * q + 2 * k + e1 * m - 1;
		VECTOR(out,ostride,to0) = w_real * sum_real - w_imag * sum_imag;
		VECTOR(out,ostride,to0 + 1) = w_real * sum_imag + w_imag * sum_real;
	      }

	    }
	}
    }


  if (q % 2 == 1)
    return 0;

  {
    double tw_arg = M_PI / ((double) factor);
    double cos_tw_arg = cos (tw_arg);
    double sin_tw_arg = sin (tw_arg);

    for (k1 = 0; k1 < product_1; k1++)
      {

	double dw_real = 1.0, dw_imag = 0.0;
	double tw_real = 1.0, tw_imag = 0.0;

	for (e1 = 0; e1 < factor; e1++)
	  {
	    double w_real, w_imag, z_real, z_imag;

	    double sum_real = 0.0;

	    if (e1 > 0)
	      {
		double tmp_real = tw_real * cos_tw_arg - tw_imag * sin_tw_arg;
		double tmp_imag = tw_real * sin_tw_arg + tw_imag * cos_tw_arg;
		tw_real = tmp_real;
		tw_imag = tmp_imag;
	      }

	    w_real = tw_real;
	    w_imag = tw_imag;

	    if (e1 > 0)
	      {
		double t_real = dw_real * cos_d_theta - dw_imag * sin_d_theta;
		double t_imag = dw_real * sin_d_theta + dw_imag * cos_d_theta;
		dw_real = t_real;
		dw_imag = t_imag;
	      }

	    for (e2 = 0; e2 <= factor - e2 - 1; e2++)
	      {

		if (e2 > 0)
		  {
		    double tmp_real = dw_real * w_real - dw_imag * w_imag;
		    double tmp_imag = dw_real * w_imag + dw_imag * w_real;
		    w_real = tmp_real;
		    w_imag = tmp_imag;
		  }


		if (e2 == factor - e2 - 1)
		  {
		    const size_t from0 = factor * k1 * q + q + 2 * e2 * q - 1;
		    z_real = VECTOR(in,istride,from0);
		    z_imag = 0.0;
		    sum_real += w_real * z_real - w_imag * z_imag;
		  }
		else
		  {
		    const size_t from0 = factor * k1 * q + q + 2 * e2 * q - 1;
		    z_real = VECTOR(in,istride,from0);
		    z_imag = VECTOR(in,istride,from0 + 1);
		    sum_real += 2 * (w_real * z_real - w_imag * z_imag);
		  }

	      }

	    {
	      const size_t to0 = k1 * q + q + e1 * m - 1;
	      VECTOR(out,ostride,to0) = sum_real;
	    }
	  }
      }
  }
  return 0;
}
