#include <config.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_histogram2d.h>

void
gsl_histogram2d_reset (gsl_histogram2d * h)
{
  size_t i;
  const size_t nx = h->nx;
  const size_t ny = h->ny;

  for (i = 0; i < nx * ny; i++)
    {
      h->bin[i] = 0;
    }
}
