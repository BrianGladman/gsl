#include <math.h>
#include <gsl_complex.h>
#include <gsl_fft_halfcomplex.h>

#include <fft_halfcomplex.h>

int
gsl_fft_halfcomplex_pass_n (const double from[],
			    double to[],
			    const unsigned int factor,
			    const unsigned int product,
			    const unsigned int n,
			    const complex twiddle[])
{

  unsigned int k, k1;

  const unsigned int m = n / factor;
  const unsigned int q = n / product;
  const unsigned int product_1 = product / factor;

  unsigned int e1, e2;

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
		  unsigned int from_idx = factor * k1 * q;
		  z_real = from[from_idx];
		  z_imag = 0.0;
		  sum_real += w_real * z_real - w_imag * z_imag;
		}
	      else if (e2 == factor - e2)
		{
		  unsigned int from_idx = factor * q * k1 + 2 * e2 * q - 1;
		  z_real = from[from_idx];
		  z_imag = 0.0;
		  sum_real += w_real * z_real;
		}
	      else
		{
		  unsigned int from_idx = factor * q * k1 + 2 * e2 * q - 1;
		  z_real = from[from_idx];
		  z_imag = from[from_idx + 1];
		  sum_real += 2 * (w_real * z_real - w_imag * z_imag);
		}

	    }

	  {
	    const unsigned int to_idx = q * k1 + e1 * m;
	    to[to_idx] = sum_real;
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
		      const unsigned int from0 = factor * k1 * q + 2 * k + 2 * e2 * q - 1;
		      z_real = from[from0];
		      z_imag = from[from0 + 1];
		    }
		  else
		    {
		      const unsigned int from0 = factor * k1 * q - 2 * k + 2 * (factor - e2) * q - 1;
		      z_real = from[from0];
		      z_imag = -from[from0 + 1];
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
		  unsigned int tskip = (q + 1) / 2 - 1;
		  w_real = twiddle[k - 1 + tskip * (e1 - 1)].real;
		  w_imag = twiddle[k - 1 + tskip * (e1 - 1)].imag;
		}

	      {
		const unsigned int to0 = k1 * q + 2 * k + e1 * m - 1;
		to[to0] = w_real * sum_real - w_imag * sum_imag;
		to[to0 + 1] = w_real * sum_imag + w_imag * sum_real;
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
		    const unsigned int from0 = factor * k1 * q + q + 2 * e2 * q - 1;
		    z_real = from[from0];
		    z_imag = 0.0;
		    sum_real += w_real * z_real - w_imag * z_imag;
		  }
		else
		  {
		    const unsigned int from0 = factor * k1 * q + q + 2 * e2 * q - 1;
		    z_real = from[from0];
		    z_imag = from[from0 + 1];
		    sum_real += 2 * (w_real * z_real - w_imag * z_imag);
		  }

	      }

	    {
	      const unsigned int to0 = k1 * q + q + e1 * m - 1;
	      to[to0] = sum_real;
	    }
	  }
      }
  }
  return 0;
}
