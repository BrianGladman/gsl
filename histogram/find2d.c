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

  return 0;
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
