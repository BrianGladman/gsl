#include <config.h>
#include <gsl_errno.h>
#include <gsl_histogram2d.h>

double
gsl_histogram2d_xmax (const gsl_histogram2d * h)
{
  const int nx = h->nx;
  return h->xrange[nx];
}

double
gsl_histogram2d_xmin (const gsl_histogram2d * h)
{
  return h->xrange[0];
}

double
gsl_histogram2d_ymax (const gsl_histogram2d * h)
{
  const int ny = h->ny;
  return h->yrange[ny];
}

double
gsl_histogram2d_ymin (const gsl_histogram2d * h)
{
  return h->yrange[0];
}

size_t
gsl_histogram2d_nx (const gsl_histogram2d * h)
{
  return h->nx;
}

size_t
gsl_histogram2d_ny (const gsl_histogram2d * h)
{
  return h->ny;
}
