#include <stdlib.h>
#include <math.h>
#include <gsl_test.h>
#include <gsl_statistics.h>

#include "test.h"

/* Test program for mean.c.  JimDavies 7.96 */

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

  double * sorted ;

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
    double sd = gsl_stats_sd (groupa, na);
    double expected = 0.0337398026922845;
    gsl_test (!within_fuzz (sd, expected),
	      "gsl_stats_sd (%g observed vs %g expected)",
	      sd, expected);
  }

  {
    double sd_est = gsl_stats_est_sd (groupa, na);
    double expected = 0.0350134479659107;
    gsl_test (!within_fuzz (sd_est, expected),
	      "gsl_stats_est_sd (%g observed vs %g expected)",
	      sd_est, expected);
  }

  {
    double absdev = gsl_stats_absdev (groupa, na);
    double expected = 0.0287571428571429;
    gsl_test (!within_fuzz (absdev, expected),
	      "gsl_stats_absdev (%g observed vs %g expected)",
	      absdev, expected);
  }

  {
    double skew = gsl_stats_skew (groupa, na);
    double expected = 0.0954642051479004;
    gsl_test (!within_fuzz (skew, expected),
	      "gsl_stats_skew (%g observed vs %g expected)",
	      skew, expected);
  }

  {
    double kurt = gsl_stats_kurtosis (groupa, na);
    double expected = -1.38583851548909 ;
    gsl_test (!within_fuzz (kurt, expected),
	      "gsl_stats_kurtosis (%g observed vs %g expected)",
	      kurt, expected);
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

  {
    int max_index = gsl_stats_max_index (groupa, na);
    int expected = 4;
    gsl_test (max_index != expected,
	      "gsl_stats_max_index (%d observed vs %d expected)",
	      max_index, expected);
  }

  {
    int min_index = gsl_stats_min_index (groupa, na);
    int expected = 3;
    gsl_test (min_index != expected,
	      "gsl_stats_min_index (%d observed vs %d expected)",
	      min_index, expected);
  }

  sorted = (double *) malloc(na * sizeof(double)) ;
  memcpy(sorted, groupa, na * sizeof(double)) ;
  
  gsl_stats_sort_data(sorted, na) ;

  {
    double median = gsl_stats_median_from_sorted_data(sorted, na) ;
    double expected = 0.07505;
    gsl_test (!within_fuzz(median,expected),
	      "gsl_stats_median_from_sorted_data (even) (%g obs vs %g exp)",
	      median, expected);
  }

  {
    double median = gsl_stats_median_from_sorted_data(sorted, na - 1) ;
    double expected = 0.0728;
    gsl_test (!within_fuzz(median,expected),
	      "gsl_stats_median_from_sorted_data (odd) (%g obs vs %g exp)",
	      median, expected);
  }


  {
    double zeroth = gsl_stats_percentile_from_sorted_data(sorted, na, 0.0) ;
    double expected = 0.0242;
    gsl_test (!within_fuzz(zeroth,expected),
	      "gsl_stats_percentile_from_sorted_data (0) (%g obs vs %g exp)",
	      zeroth, expected);
  }

  {
    double top = gsl_stats_percentile_from_sorted_data(sorted, na, 1.0) ;
    double expected = 0.1331;
    gsl_test (!within_fuzz(top,expected),
	      "gsl_stats_percentile_from_sorted_data (100) (%g obs vs %g exp)",
	      top, expected);
  }

  {
    double median = gsl_stats_percentile_from_sorted_data(sorted, na, 0.5) ;
    double expected = 0.07505;
    gsl_test (!within_fuzz(median,expected),
	      "gsl_stats_percentile_from_sorted_data (50even) (%g ob vs %g ex)",
	      median, expected);
  }

  {
    double median = gsl_stats_percentile_from_sorted_data(sorted, na - 1, 0.5);
    double expected = 0.0728;
    gsl_test (!within_fuzz(median,expected),
	      "gsl_stats_percentile_from_sorted_data (50odd) (%g obs vs %g exp)",
	      median, expected);
  }

  return gsl_test_summary ();
}

