#include <math.h>
#include <gsl_complex.h>
#include <gsl_fft_real.h>

int
gsl_fft_real_pass_n (const double from[],
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

  for (k1 = 0; k1 < q; k1++)
    {
      /* compute x = W(factor) z, for z real */

      double dw_real = 1.0, dw_imag = 0.0;

      for (e1 = 0; e1 <= factor - e1; e1++)
	{
	  double sum_real = 0.0;
	  double sum_imag = 0.0;

	  double w_real = 1.0, w_imag = 0.0;

	  if (e1 > 0)
	    {
	      double tmp_real = dw_real * cos_d_theta + dw_imag * sin_d_theta;
	      double tmp_imag = -dw_real * sin_d_theta + dw_imag * cos_d_theta;
	      dw_real = tmp_real;
	      dw_imag = tmp_imag;
	    }

	  for (e2 = 0; e2 < factor; e2++)
	    {
	      double z_real = from[k1 * product_1 + e2 * m];

	      if (e2 > 0)
		{
		  double tmp_real = dw_real * w_real - dw_imag * w_imag;
		  double tmp_imag = dw_real * w_imag + dw_imag * w_real;
		  w_real = tmp_real;
		  w_imag = tmp_imag;
		}

	      sum_real += w_real * z_real;
	      sum_imag += w_imag * z_real;

	    }
	  if (e1 == 0)
	    {
	      const unsigned int to0 = product * k1;
	      to[to0] = sum_real;
	    }
	  else if (e1 < factor - e1)
	    {
	      const unsigned int to0 = k1 * product + 2 * e1 * product_1 - 1;
	      to[to0] = sum_real;
	      to[to0 + 1] = sum_imag;
	    }
	  else if (e1 == factor - e1)
	    {
	      const unsigned int to0 = k1 * product + 2 * e1 * product_1 - 1;
	      to[to0] = sum_real;
	    }

	}
    }

  if (product_1 == 1)
    return 0;

  for (k = 1; k < (product_1 + 1) / 2; k++)
    {
      for (k1 = 0; k1 < q; k1++)
	{

	  double dw_real = 1.0, dw_imag = 0.0;

	  for (e1 = 0; e1 < factor; e1++)
	    {
	      double sum_real = 0.0, sum_imag = 0.0;

	      double w_real = 1.0, w_imag = 0.0;

	      if (e1 > 0)
		{
		  const double tmp_real = dw_real * cos_d_theta + dw_imag * sin_d_theta;
		  const double tmp_imag = -dw_real * sin_d_theta + dw_imag * cos_d_theta;
		  dw_real = tmp_real;
		  dw_imag = tmp_imag;
		}

	      for (e2 = 0; e2 < factor; e2++)
		{

		  int tskip = (product_1 + 1) / 2 - 1;
		  const unsigned int from0 = k1 * product_1 + 2 * k + e2 * m - 1;
		  double tw_real, tw_imag;
		  double z_real, z_imag;

		  if (e2 == 0)
		    {
		      tw_real = 1.0;
		      tw_imag = 0.0;
		    }
		  else
		    {
		      const unsigned int t_index = (k - 1) + (e2 - 1) * tskip;
		      tw_real = twiddle[t_index].real;
		      tw_imag = -twiddle[t_index].imag;
		    }

		  {
		    const double f0_real = from[from0];
		    const double f0_imag = from[from0 + 1];

		    z_real = tw_real * f0_real - tw_imag * f0_imag;
		    z_imag = tw_real * f0_imag + tw_imag * f0_real;
		  }

		  if (e2 > 0)
		    {
		      const double tmp_real = dw_real * w_real - dw_imag * w_imag;
		      const double tmp_imag = dw_real * w_imag + dw_imag * w_real;
		      w_real = tmp_real;
		      w_imag = tmp_imag;
		    }

		  sum_real += w_real * z_real - w_imag * z_imag;
		  sum_imag += w_real * z_imag + w_imag * z_real;
		}

	      if (e1 < factor - e1)
		{
		  const unsigned int to0 = k1 * product - 1 + 2 * e1 * product_1 + 2 * k;
		  to[to0] = sum_real;
		  to[to0 + 1] = sum_imag;
		}
	      else
		{
		  const unsigned int to0 = k1 * product - 1 + 2 * (factor - e1) * product_1 - 2 * k;
		  to[to0] = sum_real;
		  to[to0 + 1] = -sum_imag;
		}

	    }
	}
    }


  if (product_1 % 2 == 1)
    return 0;

  {
    double tw_arg = M_PI / ((double) factor);
    double cos_tw_arg = cos (tw_arg);
    double sin_tw_arg = -sin (tw_arg);

    for (k1 = 0; k1 < q; k1++)
      {
	double dw_real = 1.0, dw_imag = 0.0;

	for (e1 = 0; e1 < factor; e1++)
	  {
	    double z_real, z_imag;

	    double sum_real = 0.0;
	    double sum_imag = 0.0;

	    double w_real = 1.0, w_imag = 0.0;
	    double tw_real = 1.0, tw_imag = 0.0;

	    if (e1 > 0)
	      {
		double t_real = dw_real * cos_d_theta + dw_imag * sin_d_theta;
		double t_imag = -dw_real * sin_d_theta + dw_imag * cos_d_theta;
		dw_real = t_real;
		dw_imag = t_imag;
	      }

	    for (e2 = 0; e2 < factor; e2++)
	      {

		if (e2 > 0)
		  {
		    double tmp_real = tw_real * cos_tw_arg - tw_imag * sin_tw_arg;
		    double tmp_imag = tw_real * sin_tw_arg + tw_imag * cos_tw_arg;
		    tw_real = tmp_real;
		    tw_imag = tmp_imag;
		  }

		if (e2 > 0)
		  {
		    double tmp_real = dw_real * w_real - dw_imag * w_imag;
		    double tmp_imag = dw_real * w_imag + dw_imag * w_real;
		    w_real = tmp_real;
		    w_imag = tmp_imag;
		  }


		{
		  const unsigned int from0 = k1 * product_1 + 2 * k + e2 * m - 1;
		  const double f0_real = from[from0];
		  z_real = tw_real * f0_real;
		  z_imag = tw_imag * f0_real;
		}

		sum_real += w_real * z_real - w_imag * z_imag;
		sum_imag += w_real * z_imag + w_imag * z_real;
	      }

	    if (e1 + 1 < factor - e1)
	      {
		const unsigned int to0 = k1 * product - 1 + 2 * e1 * product_1 + 2 * k;
		to[to0] = sum_real;
		to[to0 + 1] = sum_imag;
	      }
	    else if (e1 + 1 == factor - e1)
	      {
		const unsigned int to0 = k1 * product - 1 + 2 * e1 * product_1 + 2 * k;
		to[to0] = sum_real;
	      }
	    else
	      {
		const unsigned int to0 = k1 * product - 1 + 2 * (factor - e1) * product_1 - 2 * k;
		to[to0] = sum_real;
		to[to0 + 1] = -sum_imag;
	      }

	  }
      }
  }
  return 0;
}
