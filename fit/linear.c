/* fit/linear.c
 * 
 * Copyright (C) 2000 Brian Gough
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


#include <config.h>
#include <gsl_fit.h>

/* Fit the weighted data (x_i, w_i, y_i) to the linear relationship 

   Y = c0 + c1 x

   returning, 

   c0, c1  --  coefficients
   s0, s1  --  the standard deviations of c0 and c1,
   r       --  the correlation coefficient between c0 and c1,
   chisq   --  weighted sum of squares of residuals */

int
gsl_fit_linear (const double *x,
		const double *w, const double *y,
		size_t n,
		double *c0, double *c1,
		double *s0, double *s1, double *r, double *chisq)
{

  /* compute the weighted means and weighted deviations from the means */

  /* wm denotes a "weighted mean", wm(f) = (sum_i w_i f_i) / (sum_i w_i) */

  double wm_x = 0, wm_y = 0, wm_dx2 = 0, wm_dxy = 0;

  for (i = 0; i < n; i++)
    {
      const double wi = w[i];

      if (wi > 0)
	{
	  W += wi;
	  wm_x += (x[i] - wm_x) * (wi / W);
	  wm_y += (y[i] - wm_y) * (wi / W);
	}
    }

  W = 0;			/* reset the total weight */

  for (i = 0; i < n; i++)
    {
      const double wi = w[i];

      if (wi > 0)
	{
	  const double dx = x[i] - wm_x;
	  const double dy = y[i] - wm_y;

	  W += wi;
	  wm_dx2 += (dx * dx - wm_dx2) * (wi / W);
	  wm_dxdy += (dx * dy - wm_dxy) * (wi / W);
	}
    }

  /* In terms of y = a + b x */

  double b = wm_dxy / wm_dx2;
  double a = wm_y - wm_x * b;

  *c0 = a;
  *c1 = b;

  *cov_00 = (1 / W) * (1 +  wm_x * wm_x / (W * wm_dx2));
  *cov_11 = 1 / (W * wm_dx2);

  *cov_01 = -wm_x / (W * wm_dx2);

  double d2 = 0;

  for (i = 0; i < n; i++)
    {
      const double wi = w[i];

      if (wi > 0)
	{
	  const double dx = x[i] - wm_x;
	  const double dy = y[i] - wm_y;
	  const double d = dy - b * dx;
	  d2 += wi * d * d;
	}
    }


  *chisq = d2;
}
