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

  return 0;
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

  return 0;
}
