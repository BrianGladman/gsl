void FUNCTION (test, func) (void);

void
FUNCTION (test, func) (void)
{
  /* sample sets of doubles */

  const unsigned int na = 14, nb = 14;

  const BASE groupa[] =
  {.0421, .0941, .1064, .0242, .1331,
   .0773, .0243, .0815, .1186, .0356,
   .0728, .0999, .0614, .0479};

  const BASE groupb[] =
  {.1081, .0986, .1566, .1961, .1125,
   .1942, .1079, .1021, .1583, .1673,
   .1675, .1856, .1688, .1512};

  BASE * sorted ;

#ifdef BASE_FLOAT
  double rel = 1e-6;
#else
  double rel = 1e-10;
#endif

  {
    double mean = FUNCTION(gsl_stats,mean) (groupa, na);
    double expected = 0.0728;
    gsl_test_rel (mean, expected, rel, NAME(gsl_stats) "_mean");
  }

  {
    double var = FUNCTION(gsl_stats,variance) (groupa, na);
    double expected = 0.00113837428571429;
    gsl_test_rel (var, expected, rel, NAME(gsl_stats) "_variance");
  }

  {
    double var_est = FUNCTION(gsl_stats,est_variance) (groupb, nb);
    double expected = 0.00124956615384615;
    gsl_test_rel (var_est, expected, rel, NAME(gsl_stats) "_est_variance");
  }

  {
    double sd = FUNCTION(gsl_stats,sd) (groupa, na);
    double expected = 0.0337398026922845;
    gsl_test_rel (sd, expected, rel, NAME(gsl_stats) "_sd");
  }

  {
    double sd_est = FUNCTION(gsl_stats,est_sd) (groupa, na);
    double expected = 0.0350134479659107;
    gsl_test_rel (sd_est, expected, rel, NAME(gsl_stats) "_est_sd");
  }

  {
    double absdev = FUNCTION(gsl_stats,absdev) (groupa, na);
    double expected = 0.0287571428571429;
    gsl_test_rel (absdev, expected, rel, NAME(gsl_stats) "_absdev");
  }

  {
    double skew = FUNCTION(gsl_stats,skew) (groupa, na);
    double expected = 0.0954642051479004;
    gsl_test_rel (skew, expected, rel, NAME(gsl_stats) "_skew");
  }

  {
    double kurt = FUNCTION(gsl_stats,kurtosis) (groupa, na);
    double expected = -1.38583851548909 ;
    gsl_test_rel (kurt, expected, rel, NAME(gsl_stats) "_kurtosis");
  }

  {
    double pv = FUNCTION(gsl_stats,pvariance) (groupa, groupb, na, nb);
    double expected = 0.00123775384615385;
    gsl_test_rel (pv, expected, rel, NAME(gsl_stats) "_pvariance");
  }

  {
    double t = FUNCTION(gsl_stats,ttest) (groupa, groupb, na, nb);
    double expected = -5.67026326985851;
    gsl_test_rel (t, expected, rel, NAME(gsl_stats) "_ttest");
  }

  {
    BASE max = FUNCTION(gsl_stats,max) (groupa, na);
    BASE expected = 0.1331;
    gsl_test  (max != expected,
	       NAME(gsl_stats) "_max (%g observed vs %g expected)", 
	       max, expected);
  }

  {
    BASE min = FUNCTION(gsl_stats,min) (groupa, na);
    BASE expected = 0.0242;
    gsl_test (min != expected,
	      NAME(gsl_stats) "_min (%g observed vs %g expected)", 
	      min, expected);
  }

  {
    int max_index = FUNCTION(gsl_stats,max_index) (groupa, na);
    int expected = 4;
    gsl_test (max_index != expected,
	      NAME(gsl_stats) "_max_index (%d observed vs %d expected)",
	      max_index, expected);
  }

  {
    int min_index = FUNCTION(gsl_stats,min_index) (groupa, na);
    int expected = 3;
    gsl_test (min_index != expected,
	      NAME(gsl_stats) "_min_index (%d observed vs %d expected)",
	      min_index, expected);
  }

  sorted = (BASE *) malloc(na * sizeof(BASE)) ;
  memcpy(sorted, groupa, na * sizeof(BASE)) ;
  
  FUNCTION(gsl_stats,sort_data)(sorted, na) ;

  {
    double median = FUNCTION(gsl_stats,median_from_sorted_data)(sorted, na) ;
    double expected = 0.07505;
    gsl_test_rel  (median,expected, rel,
		   NAME(gsl_stats) "_median_from_sorted_data (even)");
  }

  {
    double median = FUNCTION(gsl_stats,median_from_sorted_data)(sorted, na - 1) ;
    double expected = 0.0728;
    gsl_test_rel  (median,expected, rel,
		   NAME(gsl_stats) "_median_from_sorted_data");
  }


  {
    double zeroth = FUNCTION(gsl_stats,quantile_from_sorted_data)(sorted, na, 0.0) ;
    double expected = 0.0242;
    gsl_test_rel  (zeroth,expected, rel,
		   NAME(gsl_stats) "_quantile_from_sorted_data (0)");
  }

  {
    double top = FUNCTION(gsl_stats,quantile_from_sorted_data)(sorted, na, 1.0) ;
    double expected = 0.1331;
    gsl_test_rel  (top,expected, rel,
		   NAME(gsl_stats) "_quantile_from_sorted_data (100)");
  }

  {
    double median = FUNCTION(gsl_stats,quantile_from_sorted_data)(sorted, na, 0.5) ;
    double expected = 0.07505;
    gsl_test_rel  (median,expected, rel,
		   NAME(gsl_stats) "_quantile_from_sorted_data (50even)");
  }

  {
    double median = FUNCTION(gsl_stats,quantile_from_sorted_data)(sorted, na - 1, 0.5);
    double expected = 0.0728;
    gsl_test_rel  (median,expected, rel,
		   NAME(gsl_stats) "_quantile_from_sorted_data (50odd)");
  }
}
