/* histogram/get2d.c
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
#include <gsl/gsl_errno.h>
#include <gsl/gsl_histogram2d.h>

double
gsl_histogram2d_get (const gsl_histogram2d * h, const size_t i, const size_t j)
{
  const size_t nx = h->nx;
  const size_t ny = h->ny;

  if (i >= nx)
    {
      GSL_ERROR_RETURN ("index i lies outside valid range of 0 .. nx - 1",
			GSL_EDOM, 0);
    }

  if (j >= ny)
    {
      GSL_ERROR_RETURN ("index j lies outside valid range of 0 .. ny - 1",
			GSL_EDOM, 0);
    }

  return h->bin[i * ny + j];
}

int
gsl_histogram2d_get_xrange (const gsl_histogram2d * h, const size_t i,
			    double *xlower, double *xupper)
{
  const size_t nx = h->nx;

  if (i >= nx)
    {
      GSL_ERROR ("index i lies outside valid range of 0 .. nx - 1", GSL_EDOM);
    }

  *xlower = h->xrange[i];
  *xupper = h->xrange[i + 1];

  return GSL_SUCCESS;
}

int
gsl_histogram2d_get_yrange (const gsl_histogram2d * h, const size_t j,
			    double *ylower, double *yupper)
{
  const size_t ny = h->ny;

  if (j >= ny)
    {
      GSL_ERROR ("index j lies outside valid range of 0 .. ny - 1", GSL_EDOM);
    }

  *ylower = h->yrange[j];
  *yupper = h->yrange[j + 1];

  return GSL_SUCCESS;
}
