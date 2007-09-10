/* statistics/test_int.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Jim Davies, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <config.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_test.h>
#include <gsl/gsl_statistics_int.h>

#include "test.h"

int
main (void)
{
  /* sample sets of integers */
  
  const unsigned int ina = 20, inb = 20;

  const int test1[] = {1, 2, 3, 4, 5, 6} ;
  
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

  int * sorted ;


  {
    double mean = gsl_stats_int_mean (igroupa, ina);
    double expected = 17.0;
    gsl_test (!within_fuzz(mean,expected),
              "gsl_stats_int_mean (integer) (%g observed vs %g expected)",
              mean, expected);
  }

  {
    double mean = gsl_stats_int_mean (test1, 6);
    double expected = 3.5;
    gsl_test (!within_fuzz(mean,expected),
              "gsl_stats_int_mean (fractional) (%g observed vs %g expected)",
              mean, expected);
  }

  {
    double var = gsl_stats_int_variance (igroupa, ina);
    double expected = 13.7;
    gsl_test (!within_fuzz (var, expected),
              "gsl_stats_int_variance (%g observed vs %g expected)",
              var, expected);
  }

  {
    double var = gsl_stats_int_est_variance (igroupa, ina);
    double expected = 14.4210526315789;
    gsl_test (!within_fuzz (var, expected),
              "gsl_stats_int_est_variance (%g observed vs %g expected)",
              var, expected);
  }

  {
    double sd = gsl_stats_int_sd (igroupa, ina);
    double expected = 3.70135110466435;
    gsl_test (!within_fuzz (sd, expected),
              "gsl_stats_int_sd (%g observed vs %g expected)",
              sd, expected);
  }

  {
    double sd_est = gsl_stats_int_est_sd (igroupa, ina);
    double expected = 3.79750610685209;
    gsl_test (!within_fuzz (sd_est, expected),
              "gsl_stats_int_est_sd (%g observed vs %g expected)",
              sd_est, expected);
  }

  {
    double absdev = gsl_stats_int_absdev (igroupa, ina);
    double expected = 2.9;
    gsl_test (!within_fuzz (absdev, expected),
              "gsl_stats_int_absdev (%g observed vs %g expected)",
              absdev, expected);
  }

  {
    double skew = gsl_stats_int_skew (igroupa, ina);
    double expected = -0.909355923168064;
    gsl_test (!within_fuzz (skew, expected),
              "gsl_stats_int_skew (%g observed vs %g expected)",
              skew, expected);
  }

  {
    double kurt = gsl_stats_int_kurtosis (igroupa, ina);
    double expected = -0.233692524908094 ;
    gsl_test (!within_fuzz (kurt, expected),
              "gsl_stats_int_kurtosis (%g observed vs %g expected)",
              kurt, expected);
  }

  {
    double pv = gsl_stats_int_pvariance (igroupa, igroupb, ina, inb);
    double expected = 18.8421052631579;
    gsl_test (!within_fuzz (pv, expected),
              "gsl_stats_int_pvariance (%g observed vs %g expected)",
              pv, expected);
  }

  {
    double t = gsl_stats_int_ttest (igroupa, igroupb, ina, inb);
    double expected = -1.45701922702927;
    gsl_test (!within_fuzz (t, expected),
              "gsl_stats_int_ttest (%g observed vs %g expected)",
              t, expected);
  }

  {
    int max = gsl_stats_int_max (igroupa, ina);
    int expected = 22;
    gsl_test (max != expected,
              "gsl_stats_int_max (%d observed vs %d expected)", max, expected);
  }

  {
    int min = gsl_stats_int_min (igroupa, inb);
    int expected = 8;
    gsl_test (min != expected,
              "gsl_stats_int_min (%d observed vs %d expected)", min, expected);
  }

  {
    int max_index = gsl_stats_int_max_index (igroupa, ina);
    int expected = 9 ;
    gsl_test (max_index != expected,
              "gsl_stats_int_max_index (%d observed vs %d expected)",
              max_index, expected);
  }

  {
    int min_index = gsl_stats_int_min_index (igroupa, inb);
    int expected = 12 ;
    gsl_test (min_index != expected,
              "gsl_stats_int_min_index (%d observed vs %d expected)",
              min_index, expected);
  }

  sorted = (int *) malloc(ina * sizeof(int)) ;
  memcpy(sorted, igroupa, ina * sizeof(int)) ;

  gsl_stats_int_sort_data(sorted, ina) ;

  {
    double median = gsl_stats_int_median_from_sorted_data(sorted, ina) ;
    double expected = 18;
    gsl_test (!within_fuzz(median,expected),
              "gsl_stats_int_median_from_sorted_data (even) (%g observed vs %g expected)",
              median, expected);
  }

  {
    double median = gsl_stats_int_median_from_sorted_data(sorted, ina - 1) ;
    double expected = 18;
    gsl_test (!within_fuzz(median,expected),
              "gsl_stats_int_median_from_sorted_data (odd) (%g observed vs %g expected)",
              median, expected);
  }


  {
    double zeroth = gsl_stats_int_quantile_from_sorted_data(sorted, ina, 0.0) ;
    double expected = 8;
    gsl_test (!within_fuzz(zeroth,expected),
              "gsl_stats_quantile_from_sorted_data (0) (%g observed vs %g expected)",
              zeroth, expected);
  }

  {
    double top = gsl_stats_int_quantile_from_sorted_data(sorted, ina, 1.0) ;
    double expected = 22;
    gsl_test (!within_fuzz(top,expected),
              "gsl_stats_int_quantile_from_sorted_data (100) (%g obs vs %g exp)",
              top, expected);
  }

  {
    double median = gsl_stats_int_quantile_from_sorted_data(sorted, ina, 0.5) ;
    double expected = 18;
    gsl_test (!within_fuzz(median,expected),
              "gsl_stats_int_quantile_from_sorted_data (50, even) (%g obs vs %g exp)",
              median, expected);
  }

  {
    double median = gsl_stats_int_quantile_from_sorted_data(sorted, ina - 1, 0.5);
    double expected = 18;
    gsl_test (!within_fuzz(median,expected),
              "gsl_stats_int_quantile_from_sorted_data (50, odd) (%g obs vs %g exp)",
              median, expected);
  }


  
  exit (gsl_test_summary ());
}
