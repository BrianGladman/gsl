#include <stdio.h>
#include <gsl_test.h>
#include <gsl_statistics.h>

/* Test program for mean.c.  JimDavies 7.96 */

int 
main (void)
{

  /* sample sets of doubles */

  int na = 14, nb = 14 ;

  double groupa[] = { .0421, .0941, .1064, .0242, .1331,
		      .0773, .0243, .0815, .1186, .0356,
		      .0728, .0999, .0614, .0479 } ;
  
  double groupb[] = { .1081, .0986, .1566, .1961, .1125,
		      .1942, .1079, .1021, .1583, .1673,
		      .1675, .1856, .1688, .1512 } ;
  
  /* sample sets of integers */
  
  int ina = 20, inb = 20 ;

  int igroupa[] = { 17 , 18 , 16 , 18 , 12 , 
		    20 , 18 , 20 , 20 , 22 , 
		    20 , 10 ,  8 , 12 , 16 , 
		    16 , 18 , 20 , 18 , 21 } ;
  
  int igroupb[] = { 19 , 20 , 22 , 24 , 10 ,
		    25 , 20 , 22 , 21 , 23 ,
		    20 , 10 , 12 , 14 , 12 ,
		    20 , 22 , 24 , 23 , 17 } ;

  double mean, var, var_est, sd, sd_est, pv, t, max, min ;
  int maximum, minumum; 

  mean = gsl_stats_mean(groupa, na);
  gsl_test(mean < 0.072 || mean > 0.073, "gsl_stats_mean (%f vs XX)", mean);
  
  var = gsl_stats_variance(groupa, na);
  gsl_test(var < .001139 || var > .001135, "gsl_stats_variance (%f)", var);
  
  var_est = gsl_stats_est_variance(groupb, nb);
  gsl_test(var_est < .001 || var_est > .0013, "gsl_stats_est_variance (%f vs XX)",
	   var_est);

  sd = gsl_stats_stddev(groupa, na);
  gsl_test(sd < 0.0336 || sd > 0.0338, "gsl_stats_stddev (%f vs XX)", sd);
  
  sd_est = gsl_stats_est_stddev(groupa, na);
  gsl_test(sd_est < .034  || sd_est > .036, "gsl_stats_est_stddev (%f vs XX)",
	   sd_est);

  pv = gsl_stats_pvariance(groupa, groupb, na, nb);
  gsl_test(pv < 0.00122 || pv > 0.00124,"gsl_stats_pvariance, (%f vs XX)", pv);
  
  t = gsl_stats_ttest(groupa, groupb, na, nb);
  gsl_test (t < -5.68 || t > -5.66, "gsl_stats_ttest (%f vs XX)", t);

  max = gsl_stats_max(groupa, na);
  gsl_test (max != 0.1331, "gsl_stats_max (%f vs XX)", max);

  min = gsl_stats_min(groupa, na);
  gsl_test (max != 0.0242, "gsl_stats_min (%f vs XX)", min);


  /* integer tests */

  mean = gsl_stats_imean(igroupa, ina);
  gsl_test(mean != 17, "gsl_stats_imean (%f vs 17)", mean);

  var = gsl_stats_ivariance(igroupa, ina);
  gsl_test(var < 13.6 || var > 13.8,"gsl_stats_ivariance (%f vs XX)", var);
  
  var = gsl_stats_iest_variance(igroupa, ina);
  gsl_test(var < 14.42 || var > 14.425, "gsl_stats_iest_variance (%f vs XX)", 
	   var);
  
  sd = gsl_stats_istddev (igroupa, ina);
  gsl_test(sd < 3.6 || sd > 3.78, "gsl_stats_istddev (%f vs XX)", sd);
  
  sd_est = gsl_stats_iest_stddev(igroupa, ina);
  gsl_test(sd_est < 3.79 || sd_est > 3.8, "gsl_stats_iest_stddev (%f vs XX)",
	   sd_est);
  
  pv = gsl_stats_ipvariance(igroupa, igroupb, ina, inb);
  gsl_test(pv < 18.84 || pv > 18.85,"gsl_stats_ipvariance (%f vs XX)", pv);

  t = gsl_stats_ittest(igroupa, igroupb, ina, inb);
  gsl_test(t < -1.47 || t > -1.45, "gsl_stats_ittest (%f vs XX)", t);
  
  maximum = gsl_stats_imax(igroupa, ina);
  gsl_test(maximum != 22, "gsl_stats_imax (%d vs 22)", maximum);
  
  minumum = gsl_stats_imin(igroupa, inb);
  gsl_test(minumum != 8, "gsl_stats_imin (%d vs 8)", minumum);
  
  return 0;
}  






