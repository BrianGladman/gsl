#include <config.h>
#include <stdlib.h>
#include <math.h>

#include <gsl_test.h>
#include <gsl_statistics.h>

int
main (void)
{
  int i ;

  const int nacc1 = 3 ;
  const double numacc1[3] = { 10000001, 10000003, 10000002 } ;

  const int nacc2 = 1001 ;
  double numacc2[1001] ;

  const int nacc3 = 1001 ;
  double numacc3[1001] ;

  const int nacc4 = 1001 ;
  double numacc4[1001] ;

  numacc2[0] = 1.2 ;
  numacc3[0] = 1000000.2 ; 
  numacc4[0] = 10000000.2 ; 
 
  for (i = 1 ; i < 1000  ; i += 2) 
    {
      numacc2[i] = 1.1 ;
      numacc2[i+1] = 1.3 ;
      numacc3[i] = 1000000.1 ;
      numacc3[i+1] = 1000000.3 ;
      numacc4[i] = 10000000.1 ;
      numacc4[i+1] = 10000000.3 ;
    }
    
  {
    double mean = gsl_stats_mean (numacc1, nacc1);
    double estvar = gsl_stats_est_variance (numacc1, nacc1);
    double lag1 = gsl_stats_lag1_autocorrelation (numacc1, nacc1);

    double expected_mean = 10000002;
    double expected_estvar = 1;
    double expected_lag1 = -0.5;

    gsl_test_rel (mean, expected_mean, 1e-15, 
		  "nist-numacc1 gsl_stats_mean") ;
    gsl_test_rel (estvar, expected_estvar, 1e-15, 
		  "nist-numacc1 gsl_stats_est_variance") ;
    gsl_test_rel (lag1, expected_lag1, 1e-15, 
		  "nist-numacc1 gsl_stats_lag1_autocorrelation") ;
  }

  {
    double mean = gsl_stats_mean (numacc2, nacc2);
    double estsd = gsl_stats_est_sd (numacc2, nacc2);
    double lag1 = gsl_stats_lag1_autocorrelation (numacc2, nacc2);

    double expected_mean = 1.2;
    double expected_estsd = 0.1;
    double expected_lag1 = -0.999;

    gsl_test_rel (mean, expected_mean, 1e-15, 
		  "nist-numacc2 gsl_stats_mean") ;
    gsl_test_rel (estsd, expected_estsd, 1e-15, 
		  "nist-numacc2 gsl_stats_est_sd") ;
    gsl_test_rel (lag1, expected_lag1, 1e-15, 
		  "nist-numacc2 gsl_stats_lag1_autocorrelation") ;
  }

  {
    double mean = gsl_stats_mean (numacc3, nacc3);
    double estsd = gsl_stats_est_sd (numacc3, nacc3);
    double lag1 = gsl_stats_lag1_autocorrelation (numacc3, nacc3);

    double expected_mean = 1000000.2;
    double expected_estsd = 0.1;
    double expected_lag1 = -0.999;

    gsl_test_rel (mean, expected_mean, 1e-15, 
		  "nist-numacc3 gsl_stats_mean") ;
    gsl_test_rel (estsd, expected_estsd, 1e-9, 
		  "nist-numacc3 gsl_stats_est_sd") ;
    gsl_test_rel (lag1, expected_lag1, 1e-15, 
		  "nist-numacc3 gsl_stats_lag1_autocorrelation") ;
  }


  {
    double mean = gsl_stats_mean (numacc4, nacc4);
    double estsd = gsl_stats_est_sd (numacc4, nacc4);
    double lag1 = gsl_stats_lag1_autocorrelation (numacc4, nacc4);

    double expected_mean = 10000000.2;
    double expected_estsd = 0.1;
    double expected_lag1 = -0.999;

    gsl_test_rel (mean, expected_mean, 1e-15, 
		  "nist-numacc4 gsl_stats_mean") ;
    gsl_test_rel (estsd, expected_estsd, 1e-7, 
		  "nist-numacc4 gsl_stats_est_sd") ;
    gsl_test_rel (lag1, expected_lag1, 1e-15, 
		  "nist-numacc4 gsl_stats_lag1_autocorrelation") ;
  }


  return gsl_test_summary ();
}







