/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <stdlib.h>

#include <gsl_interp.h>

size_t
gsl_interp_bsearch (const double x_array[], double x,
		size_t index_lo,
		size_t index_hi
)
{
  size_t ilo = index_lo;
  size_t ihi = index_hi;
  while (ihi > ilo + 1)
    {
      size_t i = (ihi + ilo) / 2;
      if (x_array[i] > x)
	ihi = i;
      else
	ilo = i;
    }

  return ilo;
}
