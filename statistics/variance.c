#include <math.h>
#include <gsl_statistics.h>

double sum_of_squares (const double data[], unsigned int n) ;

double gsl_stats_variance (const double data[], unsigned int n)
{
  double sum = sum_of_squares(data, n) ;
  double variance = sum / n ; 

  return variance ; 
}

double gsl_stats_est_variance (const double data[], unsigned int n)
{
  double sum = sum_of_squares(data, n) ;
  double variance = sum / (n - 1) ; 

  return variance ; 

}

double gsl_stats_stddev (const double data[], unsigned int n)
{
  double variance = gsl_stats_variance(data, n);
  double stddev = sqrt(variance);

  return stddev;
}

double gsl_stats_est_stddev (const double data[], unsigned int n)
{
  double est_variance = gsl_stats_est_variance(data, n);
  double est_stddev = sqrt(est_variance);

  return est_stddev;
}

double sum_of_squares (const double data[], unsigned int n)
{
  /* takes an data of doubles and finds the variance */

  double sum = 0, mean ;
  int i;
  
  mean = gsl_stats_mean(data, n);   /* find the mean */

  /* find the sum of the squares */
  for (i = 0; i < n; i++)
    {
      double delta = (data[i] - mean) ;
      sum += delta*delta ;
    }

  return sum ; 
}


double gsl_stats_ivariance (const int data[], unsigned int n)
{
  
  /* takes an array of integers and finds the variance */

  double sum, the_mean;
  double square, difference;
  int i;
  
  sum = square = difference = 0;
  
  /* find the mean */
  the_mean = gsl_stats_imean(data, n);

  /* find sum of the squares */
  for (i=0; i<n; i++){
    difference = data[i] - the_mean;
    square = difference * difference;
    sum = sum + square;
  }

  /* find variance */
  return sum / n;
}


double gsl_stats_iest_variance (const int data[], unsigned int n)
{
  /* takes an data of integers and finds the variance */

  double sum, the_mean;
  double square, difference;
  int i;
  
  sum = square = difference = 0;
  
  /* find the mean */
  the_mean = gsl_stats_imean(data, n);

  /* find sum of the squares */
  for (i=0; i<n; i++){
    difference = data[i] - the_mean;
    square = difference * difference;
    sum = sum + square;
  }

  /* find variance */
  return sum / (n-1);
}


double gsl_stats_isd (const int data[], unsigned int n)
{
  /* finds the standard deviation of an data of integers */

  double variance, sd;

  variance = gsl_stats_ivariance(data, n);

  sd = sqrt(variance);

  return sd;
}


double gsl_stats_iest_sd (const int data[], unsigned int n)
{
  /* finds the standard deviation of an data of integers */

  double variance, sd;

  variance = gsl_stats_iest_variance(data, n);

  sd = sqrt(variance);

  return sd;
}





