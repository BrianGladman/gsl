/* fft/real_pass_n.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

static void
FUNCTION(fft_real,pass_n) (const BASE in[],
			   const size_t istride,
			   BASE out[],
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
	      double z_real = VECTOR(in,istride,k1 * product_1 + e2 * m);

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
	      const size_t to0 = product * k1;
	      VECTOR(out,ostride,to0) = sum_real;
	    }
	  else if (e1 < factor - e1)
	    {
	      const size_t to0 = k1 * product + 2 * e1 * product_1 - 1;
	      VECTOR(out,ostride,to0) = sum_real;
	      VECTOR(out,ostride,to0 + 1) = sum_imag;
	    }
	  else if (e1 == factor - e1)
	    {
	      const size_t to0 = k1 * product + 2 * e1 * product_1 - 1;
	      VECTOR(out,ostride,to0) = sum_real;
	    }

	}
    }

  if (product_1 == 1)
    return;

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
		  const size_t from0 = k1 * product_1 + 2 * k + e2 * m - 1;
		  double tw_real, tw_imag;
		  double z_real, z_imag;

		  if (e2 == 0)
		    {
		      tw_real = 1.0;
		      tw_imag = 0.0;
		    }
		  else
		    {
		      const size_t t_index = (k - 1) + (e2 - 1) * tskip;
		      tw_real = GSL_REAL(twiddle[t_index]);
		      tw_imag = -GSL_IMAG(twiddle[t_index]);
		    }

		  {
		    const double f0_real = VECTOR(in,istride,from0);
		    const double f0_imag = VECTOR(in,istride,from0 + 1);

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
		  const size_t to0 = k1 * product - 1 + 2 * e1 * product_1 + 2 * k;
		  VECTOR(out,ostride,to0) = sum_real;
		  VECTOR(out,ostride,to0 + 1) = sum_imag;
		}
	      else
		{
		  const size_t to0 = k1 * product - 1 + 2 * (factor - e1) * product_1 - 2 * k;
		  VECTOR(out,ostride,to0) = sum_real;
		  VECTOR(out,ostride,to0 + 1) = -sum_imag;
		}

	    }
	}
    }


  if (product_1 % 2 == 1)
    return;

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
		  const size_t from0 = k1 * product_1 + 2 * k + e2 * m - 1;
		  const double f0_real = VECTOR(in,istride,from0);
		  z_real = tw_real * f0_real;
		  z_imag = tw_imag * f0_real;
		}

		sum_real += w_real * z_real - w_imag * z_imag;
		sum_imag += w_real * z_imag + w_imag * z_real;
	      }

	    if (e1 + 1 < factor - e1)
	      {
		const size_t to0 = k1 * product - 1 + 2 * e1 * product_1 + 2 * k;
		VECTOR(out,ostride,to0) = sum_real;
		VECTOR(out,ostride,to0 + 1) = sum_imag;
	      }
	    else if (e1 + 1 == factor - e1)
	      {
		const size_t to0 = k1 * product - 1 + 2 * e1 * product_1 + 2 * k;
		VECTOR(out,ostride,to0) = sum_real;
	      }
	    else
	      {
		const size_t to0 = k1 * product - 1 + 2 * (factor - e1) * product_1 - 2 * k;
		VECTOR(out,ostride,to0) = sum_real;
		VECTOR(out,ostride,to0 + 1) = -sum_imag;
	      }

	  }
      }
  }
  return;
}
