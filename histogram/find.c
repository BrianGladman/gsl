/* histogram/find.c
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

typedef int cmp_fn_t (const void *, const void *);
static int compare_range (const double *x, const double *range);

int
gsl_histogram_find (const gsl_histogram * h,
		    const double x, size_t * i)
{
  int status = gsl_histogram_find_impl (h->n, h->range, x, i);

  if (status)
    {
      GSL_ERROR ("x not found in range of h", GSL_EDOM);
    }

  return GSL_SUCCESS;
}

int
gsl_histogram_find_impl (const size_t n, const double range[],
			 const double x, size_t * i)
{
  if (x < range[0])
    {
      return -1;
    }

  if (x >= range[n])
    {
      return +1;
    }

  {
    double *p = (double *) bsearch (&x, range, n, sizeof (double),
				    (cmp_fn_t *) compare_range);
    if (p == 0)
      {
	GSL_ERROR ("x not found in range", GSL_ESANITY);
      }

    *i = (p - range);
  }

  return 0;
}

static int
compare_range (const double *x, const double *range)
{
  if (*x < *range)
    {
      return -1;
    }
  else if (*x >= *(range + 1))
    {
      return +1;
    }
  else
    {
      return 0;
    }
}
