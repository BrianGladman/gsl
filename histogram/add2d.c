#include <config.h>
#include <gsl_errno.h>
#include <gsl_histogram.h>
#include <gsl_histogram2d.h>

int
gsl_histogram2d_increment (gsl_histogram2d * h, double x, double y)
{
  int status = gsl_histogram2d_accumulate (h, x, y, 1.0);
  return status;
}

int
gsl_histogram2d_accumulate (gsl_histogram2d * h,
			    double x, double y, double weight)
{
  const size_t nx = h->nx;
  const size_t ny = h->ny;

  size_t i = 0, j = 0;

  int status = gsl_histogram2d_find_impl (h, x, y, &i, &j);

  if (status)
    {
      return GSL_EDOM;
    }

  if (i >= nx)
    {
      GSL_ERROR ("index lies outside valid range of 0 .. nx - 1",
		 GSL_ESANITY);
    }

  if (j >= ny)
    {
      GSL_ERROR ("index lies outside valid range of 0 .. ny - 1",
		 GSL_ESANITY);
    }

  h->bin[i * ny + j] += weight;

  return 0;
}
