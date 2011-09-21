/* bspline/greville.c
 *
 * Copyright (C) 2006, 2007, 2008, 2009 Patrick Alken
 * Copyright (C) 2008, 2011 Rhys Ulerich
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <config.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_statistics.h>

/* Return the location of the i-th Greville abscissa */
double
gsl_bspline_greville_abscissa(size_t i, gsl_bspline_workspace *w)
{
  const size_t stride = w->knots->stride;
  size_t km1 = w->km1;
  double * data = w->knots->data + (i+1)*stride;
#if GSL_RANGE_CHECK
  if (GSL_RANGE_COND(i >= gsl_bspline_ncoeffs(w)))
    {
      GSL_ERROR_VAL ("Greville abscissa index out of range", GSL_EINVAL, 0);
    }
#endif

  if (km1 == 0)
    {
      /* Return interval midpoints in degenerate k = 1 case*/
      km1   = 2;
      data -= stride;
    }

  return gsl_stats_mean(data, stride, km1);
}
