#include <math.h>
#include <gsl_test.h>
#include <gsl_statistics.h>

/* Test program for mean.c.  JimDavies 7.96 */

int within_fuzz (double x, double y);	/* approximate comparison function */

int
main (void)
{

  /* sample sets of doubles */

  const unsigned int na = 14, nb = 14;

  const double groupa[] =
  {.0421, .0941, .1064, .0242, .1331,
   .0773, .0243, .0815, .1186, .0356,
   .0728, .0999, .0614, .0479};

  const double groupb[] =
  {.1081, .0986, .1566, .1961, .1125,
   .1942, .1079, .1021, .1583, .1673,
   .1675, .1856, .1688, .1512};

  /* sample sets of integers */

  const unsigned int ina = 20, inb = 20;

  const int igroupa[] =
  {17, 18, 16, 18, 12,
   20, 18, 20, 20, 22,
   20, 10, 8, 12, 16,
   16, 18, 20, 18, 21};

  const int igroupb[] =
  {19, 20, 22, 24, 10,
   25, 20, 22, 21, 23,
   20, 10, 12, 14, 12,
   20, 22, 24, 23, 17};

  {
    double mean = gsl_stats_mean (groupa, na);
    double expected = 0.0728;
    gsl_test (!within_fuzz (mean, expected),
	      "gsl_stats_mean (%g observed vs %g expected)",
	      mean, expected);
  }

  {
    double var = gsl_stats_variance (groupa, na);
    double expected = 0.00113837428571429;
    gsl_test (!within_fuzz (var, expected),
	      "gsl_stats_variance (%g observed vs %g expected)",
	      var, expected);
  }

  {
    double var_est = gsl_stats_est_variance (groupb, nb);
    double expected = 0.00124956615384615;
    gsl_test (!within_fuzz (var_est, expected),
	      "gsl_stats_est_variance (%g observed vs %g expected)",
	      var_est, expected);
  }

  {
    double sd = gsl_stats_stddev (groupa, na);
    double expected = 0.0337398026922845;
    gsl_test (!within_fuzz (sd, expected),
	      "gsl_stats_stddev (%g observed vs %g expected)",
	      sd, expected);
  }

  {
    double sd_est = gsl_stats_est_stddev (groupa, na);
    double expected = 0.0350134479659107;
    gsl_test (!within_fuzz (sd_est, expected),
	      "gsl_stats_est_stddev (%g observed vs %g expected)",
	      sd_est, expected);
  }

  {
    double pv = gsl_stats_pvariance (groupa, groupb, na, nb);
    double expected = 0.00123775384615385;
    gsl_test (!within_fuzz (pv, expected),
	      "gsl_stats_pvariance (%g observed vs %g expected)",
	      pv, expected);
  }

  {
    double t = gsl_stats_ttest (groupa, groupb, na, nb);
    double expected = -5.67026326985851;
    gsl_test (!within_fuzz (t, expected),
	      "gsl_stats_ttest (%g observed vs %g expected)",
	      t, expected);
  }

  {
    double max = gsl_stats_max (groupa, na);
    double expected = 0.1331;
    gsl_test (max != expected,
	      "gsl_stats_max (%g observed vs %g expected)", max, expected);
  }

  {
    double min = gsl_stats_min (groupa, na);
    double expected = 0.0242;
    gsl_test (min != expected,
	      "gsl_stats_min (%g observed vs %g expected)", min, expected);
  }

  /* integer tests */

  {
    double mean = gsl_stats_imean (igroupa, ina);
    double expected = 17;
    gsl_test (mean != expected,
	      "gsl_stats_imean (%g observed vs %g expected)",
	      mean, expected);
  }

  {
    double var = gsl_stats_ivariance (igroupa, ina);
    double expected = 13.7;
    gsl_test (!within_fuzz (var, expected),
	      "gsl_stats_ivariance (%g observed vs %g expected)",
	      var, expected);
  }

  {
    double var = gsl_stats_iest_variance (igroupa, ina);
    double expected = 14.4210526315789;
    gsl_test (!within_fuzz (var, expected),
	      "gsl_stats_iest_variance (%g observed vs %g expected)",
	      var, expected);
  }

  {
    double sd = gsl_stats_istddev (igroupa, ina);
    double expected = 3.70135110466435;
    gsl_test (!within_fuzz (sd, expected),
	      "gsl_stats_istddev (%g observed vs %g expected)",
	      sd, expected);
  }

  {
    double sd_est = gsl_stats_iest_stddev (igroupa, ina);
    double expected = 3.79750610685209;
    gsl_test (!within_fuzz (sd_est, expected),
	      "gsl_stats_iest_stddev (%g observed vs %g expected)",
	      sd_est, expected);
  }

  {
    double pv = gsl_stats_ipvariance (igroupa, igroupb, ina, inb);
    double expected = 18.8421052631579;
    gsl_test (!within_fuzz (pv, expected),
	      "gsl_stats_ipvariance (%g observed vs %g expected)",
	      pv, expected);
  }

  {
    double t = gsl_stats_ittest (igroupa, igroupb, ina, inb);
    double expected = -1.45701922702927;
    gsl_test (!within_fuzz (t, expected),
	      "gsl_stats_ittest (%g observed vs %g expected)",
	      t, expected);
  }

  {
    int max = gsl_stats_imax (igroupa, ina);
    int expected = 22;
    gsl_test (max != expected,
	      "gsl_stats_imax (%d observed vs %d expected)", max, expected);
  }

  {
    int min = gsl_stats_imin (igroupa, inb);
    int expected = 8;
    gsl_test (min != expected,
	      "gsl_stats_imin (%d observed vs %d expected)", min, expected);
  }

  return gsl_test_summary ();
}


int 
within_fuzz (double x, double y)
{
  return fabs (x - y) < 0.00001;
}
