/* histogram/find2d.c
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

#include <config.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_histogram2d.h>

int
gsl_histogram2d_find (const gsl_histogram2d * h,
		      const double x, const double y,
		      size_t * i, size_t * j)
{
  int status = gsl_histogram_find_impl (h->nx, h->xrange, x, i);

  if (status)
    {
      GSL_ERROR ("x not found in range of h", GSL_EDOM);
    }

  status = gsl_histogram_find_impl (h->ny, h->yrange, y, j);

  if (status)
    {
      GSL_ERROR ("y not found in range of h", GSL_EDOM);
    }

  return GSL_SUCCESS;
}


int
gsl_histogram2d_find_impl (const gsl_histogram2d * h,
			   const double x, const double y,
			   size_t * i, size_t * j)
{
  int status = gsl_histogram_find_impl (h->nx, h->xrange, x, i);

  if (status)
    {
      return status;
    }

  status = gsl_histogram_find_impl (h->ny, h->yrange, y, j);

  if (status)
    {
      return status;
    }

  return 0;
}
