/* multimin/diff.c
 * 
 * Copyright (C) 2000 David Morrison
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
#include <gsl/gsl_multimin.h>

int
gsl_multimin_diff (const gsl_multimin_function * f,
		   const gsl_vector * x, gsl_vector * g)
{
  /* This is nearly identical to gsl_diff_central except that the
     multimin functions are used to restrict a function requiring
     vecor input to a single real valued input. */
  size_t i, k, n;
  gsl_vector *direction;
  gsl_multimin_to_single *w;
  double h = GSL_SQRT_DBL_EPSILON;
  double a[4], d[4], a3;
  double fl, fh;

  direction = gsl_vector_calloc (f->n);
  w = gsl_multimin_to_single_alloc (f, x, direction);

  for (n = 0; n < f->n; n++)
    {
      gsl_vector_set_basis (direction, n);

      /* Algorithm based on description on pg. 204 of Conte and de Boor
         (CdB) - coefficients of Newton form of polynomial of degree 3. */

      for (i = 0; i < 4; i++)
	{
	  a[i] = (i - 2.0) * h;
	  gsl_multimin_compute_ep (w, a[i]);
	  d[i] = f->f (w->evaluation_point, 0);
	  /*      d[i] = gsl_multimin_to_single_eval(a[i], w); */
	}

      for (k = 1; k < 5; k++)
	{
	  for (i = 0; i < 4 - k; i++)
	    {
	      d[i] = (d[i + 1] - d[i]) / (a[i + k] - a[i]);
	    }
	}

      /* Adapt procedure described on pg. 282 of CdB to find best
         value of step size. */

      a3 = fabs (d[0] + d[1] + d[2] + d[3]);

      if (a3 < 100.0 * GSL_SQRT_DBL_EPSILON)
	{
	  a3 = 100.0 * GSL_SQRT_DBL_EPSILON;
	}

      h = pow (GSL_SQRT_DBL_EPSILON / (2.0 * a3), 1.0 / 3.0);

      if (h > 100.0 * GSL_SQRT_DBL_EPSILON)
	{
	  h = 100.0 * GSL_SQRT_DBL_EPSILON;
	}

      gsl_multimin_compute_ep (w, h);
      fh = f->f (w->evaluation_point, 0);
      gsl_multimin_compute_ep (w, -h);
      fl = f->f (w->evaluation_point, 0);
      gsl_vector_set (g, n, (fh - fl) / (2.0 * h));
    }

  gsl_vector_free (direction);
  gsl_multimin_to_single_free (w);

  return GSL_SUCCESS;
}
