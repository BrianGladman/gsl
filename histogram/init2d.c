#include <config.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_histogram2d.h>

gsl_histogram2d *
gsl_histogram2d_calloc_uniform (const size_t nx, const size_t ny,
				const double xmin, const double xmax,
				const double ymin, const double ymax)
{
  gsl_histogram2d *h;

  if (xmin >= xmax)
    {
      GSL_ERROR_RETURN ("xmin must be less than xmax", GSL_EINVAL, 0);
    }

  if (ymin >= ymax)
    {
      GSL_ERROR_RETURN ("ymin must be less than ymax", GSL_EINVAL, 0);
    }

  h = gsl_histogram2d_calloc (nx, ny);

  if (h == 0)
    {
      return h;
    }

  {
    size_t i;

    for (i = 0; i < nx + 1; i++)
      {
	h->xrange[i] = xmin + ((double) i / (double) nx) * (xmax - xmin);
      }

    for (i = 0; i < ny + 1; i++)
      {
	h->yrange[i] = ymin + ((double) i / (double) ny) * (ymax - ymin);
      }
  }

  return h;
}

gsl_histogram2d *
gsl_histogram2d_calloc (const size_t nx, const size_t ny)
{
  gsl_histogram2d *h;

  if (nx == 0)
    {
      GSL_ERROR_RETURN ("histogram2d length nx must be positive integer",
			GSL_EDOM, 0);
    }

  if (ny == 0)
    {
      GSL_ERROR_RETURN ("histogram2d length ny must be positive integer",
			GSL_EDOM, 0);
    }

  h = (gsl_histogram2d *) malloc (sizeof (gsl_histogram2d));

  if (h == 0)
    {
      GSL_ERROR_RETURN ("failed to allocate space for histogram2d struct",
			GSL_ENOMEM, 0);
    }

  h->xrange = (double *) malloc ((nx + 1) * sizeof (double));

  if (h->xrange == 0)
    {
      free (h);		/* exception in constructor, avoid memory leak */

      GSL_ERROR_RETURN ("failed to allocate space for histogram2d x ranges",
			GSL_ENOMEM, 0);
    }

  h->yrange = (double *) malloc ((ny + 1) * sizeof (double));

  if (h->yrange == 0)
    {
      free (h->xrange);
      free (h);		/* exception in constructor, avoid memory leak */

      GSL_ERROR_RETURN ("failed to allocate space for histogram2d y ranges",
			GSL_ENOMEM, 0);
    }

  h->bin = (double *) malloc (nx * ny * sizeof (double));

  if (h->bin == 0)
    {
      free (h->xrange);
      free (h->yrange);
      free (h);		/* exception in constructor, avoid memory leak */

      GSL_ERROR_RETURN ("failed to allocate space for histogram bins",
			GSL_ENOMEM, 0);
    }

  {
    size_t i;

    for (i = 0; i < nx + 1; i++)
      {
	h->xrange[i] = i;
      }

    for (i = 0; i < ny + 1; i++)
      {
	h->yrange[i] = i;
      }

    for (i = 0; i < nx * ny; i++)
      {
	h->bin[i] = 0;
      }
  }

  h->nx = nx;
  h->ny = ny;

  return h;
}


void
gsl_histogram2d_free (gsl_histogram2d * h)
{
  free (h->xrange);
  free (h->yrange);
  free (h->bin);
  free (h);
}
