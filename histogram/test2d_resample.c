#include <config.h>
#include <math.h>
#include <gsl_histogram2d.h>
#include <gsl_test.h>

#include "urand.c"

int
main (void)
{
  size_t i, j;
  int status = 0;
  double total = 0;
  size_t N = 200000;

  gsl_histogram2d *h = gsl_histogram2d_calloc_uniform (10, 10,
						       0.0, 1.0,
						       0.0, 1.0);
  for (i = 0; i < 10; i++)
    {
      for (j = 0; j < 10; j++)
	{
	  double w = 10.0 * i + j;
	  total += w;
	  gsl_histogram2d_accumulate (h, 0.1 * i, 0.1 * i, w);
	}
    }

  {
    gsl_histogram2d_pdf *p = gsl_histogram2d_pdf_alloc (h);

    gsl_histogram2d *hh = gsl_histogram2d_calloc_uniform (20, 20,
							  0.0, 1.0,
							  0.0, 1.0);

    for (i = 0; i < N; i++)
      {
	double u = urand();
	double v = urand();
	double x, y;
	status = gsl_histogram2d_pdf_sample (p, u, v, &x, &y);
	status = gsl_histogram2d_increment (hh, x, y);
      }

    status = 0;
    for (i = 0; i < 20; i++)
      {
	for (j = 0; j < 20; j++)
	  {
	    double z = 4 * total * gsl_histogram2d_get (hh, i, j) / (double) N;
	    size_t k1, k2;
	    double ya;
	    double x, xmax, y, ymax;

	    gsl_histogram2d_get_xrange (hh, i, &x, &xmax);
	    gsl_histogram2d_get_yrange (hh, j, &y, &ymax);

	    gsl_histogram2d_find (h, x, y, &k1, &k2);
	    ya = gsl_histogram2d_get (h, k1, k2);

	    if (ya == 0)
	      {
		if (z != 0)
		  {
		    status = 1;
		    printf ("(%d,%d): %g vs %g\n", (int)i, (int)j, z, ya);
		  }
	      }
	    else
	      {
		double err = 1 / sqrt (gsl_histogram2d_get (hh, i, j));
		double sigma = fabs ((z - ya) / (ya * err));
		if (sigma > 3)
		  {
		    status = 1;
		    printf ("%g vs %g err=%g sigma=%g\n", z, ya, err, sigma);
		  }
	      }
	  }
      }
    gsl_test (status, "gsl_histogram2d_pdf_sample within statistical errors");
  }
  return gsl_test_summary ();
}
