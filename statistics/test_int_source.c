void FUNCTION (test, func) (void);

void
FUNCTION (test, func) (void)
{
  /* sample sets of integers */
  
  const unsigned int ina = 20, inb = 20;

  const BASE test1[] = {1, 2, 3, 4, 5, 6} ;
  
  const BASE igroupa[] =
  {17, 18, 16, 18, 12,
   20, 18, 20, 20, 22,
   20, 10, 8, 12, 16,
   16, 18, 20, 18, 21};

  const BASE igroupb[] =
  {19, 20, 22, 24, 10,
   25, 20, 22, 21, 23,
   20, 10, 12, 14, 12,
   20, 22, 24, 23, 17};

  BASE * sorted ;

  double rel = 1e-10 ;

  {
    double mean = FUNCTION(gsl_stats,mean) (igroupa, ina);
    double expected = 17.0;
    gsl_test_rel (mean,expected, rel, NAME(gsl_stats) "_mean (integer)");
  }

  {
    double mean = FUNCTION(gsl_stats,mean) (test1, 6);
    double expected = 3.5;
    gsl_test_rel (mean,expected, rel, NAME(gsl_stats) "_mean (fractional)");
  }

  {
    double var = FUNCTION(gsl_stats,variance) (igroupa, ina);
    double expected = 13.7;
    gsl_test_rel (var, expected, rel, NAME(gsl_stats) "_variance");
  }

  {
    double var = FUNCTION(gsl_stats,est_variance) (igroupa, ina);
    double expected = 14.4210526315789;
    gsl_test_rel (var, expected, rel, NAME(gsl_stats) "_est_variance");
  }

  {
    double sd = FUNCTION(gsl_stats,sd) (igroupa, ina);
    double expected = 3.70135110466435;
    gsl_test_rel (sd, expected, rel, NAME(gsl_stats) "_sd");
  }

  {
    double sd_est = FUNCTION(gsl_stats,est_sd) (igroupa, ina);
    double expected = 3.79750610685209;
    gsl_test_rel (sd_est, expected, rel, NAME(gsl_stats) "_est_sd");
  }

  {
    double absdev = FUNCTION(gsl_stats,absdev) (igroupa, ina);
    double expected = 2.9;
    gsl_test_rel (absdev, expected, rel, NAME(gsl_stats) "_absdev");
  }

  {
    double skew = FUNCTION(gsl_stats,skew) (igroupa, ina);
    double expected = -0.909355923168064;
    gsl_test_rel (skew, expected, rel, NAME(gsl_stats) "_skew");
  }

  {
    double kurt = FUNCTION(gsl_stats,kurtosis) (igroupa, ina);
    double expected = -0.233692524908094 ;
    gsl_test_rel (kurt, expected, rel, NAME(gsl_stats) "_kurtosis");
  }

  {
    double pv = FUNCTION(gsl_stats,pvariance) (igroupa, igroupb, ina, inb);
    double expected = 18.8421052631579;
    gsl_test_rel (pv, expected, rel, NAME(gsl_stats) "_pvariance");
  }

  {
    double t = FUNCTION(gsl_stats,ttest) (igroupa, igroupb, ina, inb);
    double expected = -1.45701922702927;
    gsl_test_rel (t, expected, rel, NAME(gsl_stats) "_ttest");
  }

  {
    int max = FUNCTION(gsl_stats,max) (igroupa, ina);
    int expected = 22;
    gsl_test (max != expected,
	      NAME(gsl_stats) "_max (%d observed vs %d expected)", max, expected);
  }

  {
    int min = FUNCTION(gsl_stats,min) (igroupa, inb);
    int expected = 8;
    gsl_test (min != expected,
	      NAME(gsl_stats) "_min (%d observed vs %d expected)", min, expected);
  }

  {
    int max_index = FUNCTION(gsl_stats,max_index) (igroupa, ina);
    int expected = 9 ;
    gsl_test (max_index != expected,
	      NAME(gsl_stats) "_max_index (%d observed vs %d expected)",
	      max_index, expected);
  }

  {
    int min_index = FUNCTION(gsl_stats,min_index) (igroupa, inb);
    int expected = 12 ;
    gsl_test (min_index != expected,
	      NAME(gsl_stats) "_min_index (%d observed vs %d expected)",
	      min_index, expected);
  }

  sorted = (BASE *) malloc(ina * sizeof(BASE)) ;
  memcpy(sorted, igroupa, ina * sizeof(BASE)) ;

  FUNCTION(gsl_stats,sort_data(sorted, ina)) ;

  {
    double median = FUNCTION(gsl_stats,median_from_sorted_data)(sorted, ina) ;
    double expected = 18;
    gsl_test_rel (median,expected, rel,
		  NAME(gsl_stats) "_median_from_sorted_data (even)");
  }

  {
    double median = FUNCTION(gsl_stats,median_from_sorted_data)(sorted, ina - 1) ;
    double expected = 18;
    gsl_test_rel (median,expected, rel,
		  NAME(gsl_stats) "_median_from_sorted_data (odd)");
  }


  {
    double zeroth = FUNCTION(gsl_stats,quantile_from_sorted_data)(sorted, ina, 0.0) ;
    double expected = 8;
    gsl_test_rel (zeroth,expected, rel,
		  NAME(gsl_stats) "_quantile_from_sorted_data (0)");
  }

  {
    double top = FUNCTION(gsl_stats,quantile_from_sorted_data)(sorted, ina, 1.0) ;
    double expected = 22;
    gsl_test_rel (top,expected, rel,
		  NAME(gsl_stats) "_quantile_from_sorted_data (100)");
  }

  {
    double median = FUNCTION(gsl_stats,quantile_from_sorted_data)(sorted, ina, 0.5) ;
    double expected = 18;
    gsl_test_rel (median,expected, rel,
		  NAME(gsl_stats) "_quantile_from_sorted_data (50, even)");
  }

  {
    double median = FUNCTION(gsl_stats,quantile_from_sorted_data)(sorted, ina - 1, 0.5);
    double expected = 18;
    gsl_test_rel (median,expected, rel,
		  NAME(gsl_stats) "_quantile_from_sorted_data (50, odd)");
  }

}
