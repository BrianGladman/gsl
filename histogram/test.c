#include <config.h>
#include <stdio.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_test.h>

#define N 397

int
main (void)
{
  gsl_histogram *h;
  size_t i, j;

  h = gsl_histogram_calloc (N);

  gsl_test (h->range == 0, "gsl_histogram_calloc returns valid range pointer");
  gsl_test (h->bin == 0, "gsl_histogram_calloc returns valid bin pointer");
  gsl_test (h->n != N, "gsl_histogram_calloc returns valid size");

  for (i = 0; i < N; i++)
    {
      gsl_histogram_accumulate (h, (double) i, (double) i);
    };

  {
    int status = 0;

    for (i = 0; i < N; i++)
      {
	if (h->bin[i] != (double) i)
	  {
	    status = 1;
	  }
      };

    gsl_test (status, "gsl_histogram_accumulate writes into array correctly");
  }


  {
    int status = 0;

    for (i = 0; i < N; i++)
      {
	if (gsl_histogram_get (h, i) != i)
	  status = 1;
      };
    gsl_test (status, "gsl_histogram_get reads from array correctly");
  }

  gsl_histogram_reset (h);

  {
    int status = 0;

    for (i = 0; i < N; i++)
      {
	if (h->bin[i] != 0)
	  status = 1;
      }
    gsl_test (status, "gsl_histogram_reset zeros array correctly");
  }


  {
    int status = 0;

    for (i = 0; i < N; i++)
      {
	gsl_histogram_increment (h, (double) i);

	for (j = 0; j <= i; j++)
	  {
	    if (h->bin[j] != 1)
	      {
		status = 1;
	      }
	  }

	for (j = i + 1; j < N; j++)
	  {
	    if (h->bin[j] != 0)
	      {
		status = 1;
	      }
	  }
      }

    gsl_test (status, "gsl_histogram_increment works correctly");
  }

  {
    int status = 0;
    for (i = 0; i < N; i++)
      {
	double x0 = 0, x1 = 0;

	gsl_histogram_get_range (h, i, &x0, &x1);

	if (x0 != i || x1 != i + 1)
	  {
	    status = 1;
	  }
      }
    gsl_test (status, "gsl_histogram_getbinrange works correctly");
  }

  {
    int status = 0;
    if (gsl_histogram_max (h) != N)
      status = 1;
    gsl_test (status, "gsl_histogram_max works correctly");
  }

  {
    int status = 0;
    if (gsl_histogram_min (h) != 0)
      status = 1;
    gsl_test (status, "gsl_histogram_min works correctly");
  }

  {
    int status = 0;
    if (gsl_histogram_bins (h) != N)
      status = 1;
    gsl_test (status, "gsl_histogram_bins works correctly");
  }

  gsl_histogram_free (h);	/* free whatever is in h */

  h = gsl_histogram_calloc_uniform (N, 0.0, 1.0);

  gsl_test (h->range == 0,
	    "gsl_histogram_calloc_uniform returns valid range pointer");
  gsl_test (h->bin == 0,
	    "gsl_histogram_calloc_uniform returns valid bin pointer");
  gsl_test (h->n != N,
	    "gsl_histogram_calloc_uniform returns valid size");

  gsl_histogram_accumulate (h, 0.0, 1.0);
  gsl_histogram_accumulate (h, 0.1, 2.0);
  gsl_histogram_accumulate (h, 0.2, 3.0);
  gsl_histogram_accumulate (h, 0.3, 4.0);

  {
    size_t i1, i2, i3, i4;
    double expected;
    int status = gsl_histogram_find (h, 0.0, &i1);
    status = gsl_histogram_find (h, 0.1, &i2);
    status = gsl_histogram_find (h, 0.2, &i3);
    status = gsl_histogram_find (h, 0.3, &i4);

    for (i = 0; i < N; i++)
      {
	if (i == i1)
	  {
	    expected = 1.0;
	  }
	else if (i == i2)
	  {
	    expected = 2.0;
	  }
	else if (i == i3)
	  {
	    expected = 3.0;
	  }
	else if (i == i4)
	  {
	    expected = 4.0;
	  }
	else
	  {
	    expected = 0.0;
	  }

	if (h->bin[i] != expected)
	  {
	    status = 1;
	  }

      }
    gsl_test (status, "gsl_histogram_find works correctly");
  }


  {
    FILE *f = fopen ("test.txt", "w");
    gsl_histogram_fprintf (f, h, "%.19g", "%.19g");
    fclose (f);
  }

  {
    FILE *f = fopen ("test.txt", "r");
    gsl_histogram *hh = gsl_histogram_calloc (N);
    int status = 0;

    gsl_histogram_fscanf (f, hh);

    for (i = 0; i < N; i++)
      {
	if (h->range[i] != hh->range[i])
	  status = 1;
	if (h->bin[i] != hh->bin[i])
	  status = 1;
      }
    if (h->range[N] != hh->range[N])
      status = 1;

    gsl_test (status, "gsl_histogram_fprintf and fscanf work correctly");

    gsl_histogram_free (hh);
    fclose (f);
  }

  {
    FILE *f = fopen ("test.dat", "w");
    gsl_histogram_fwrite (f, h);
    fclose (f);
  }

  {
    FILE *f = fopen ("test.dat", "r");
    gsl_histogram *hh = gsl_histogram_calloc (N);
    int status = 0;

    gsl_histogram_fread (f, hh);

    for (i = 0; i < N; i++)
      {
	if (h->range[i] != hh->range[i])
	  status = 1;
	if (h->bin[i] != hh->bin[i])
	  status = 1;
      }
    if (h->range[N] != hh->range[N])
      status = 1;

    gsl_test (status, "gsl_histogram_fwrite and fread work correctly");

    gsl_histogram_free (hh);
    fclose (f);
  }

  gsl_histogram_free (h);

  return gsl_test_summary ();
}
