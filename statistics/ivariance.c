#include <math.h>
#include <gsl_statistics.h>

#include "ivariance.h"

double 
gsl_stats_ivariance (const int data[], unsigned int n)
{
  double mean = gsl_stats_imean (data, n);
  return gsl_stats_ivariance_with_mean(data, n, mean);
}

double 
gsl_stats_iest_variance (const int data[], unsigned int n)
{
  double mean = gsl_stats_imean (data, n);
  return gsl_stats_iest_variance_with_mean(data, n, mean);

}

double 
gsl_stats_istddev (const int data[], unsigned int n)
{
  double mean = gsl_stats_imean (data, n);
  return gsl_stats_istddev_with_mean (data, n, mean) ;
}

double 
gsl_stats_iest_stddev (const int data[], unsigned int n)
{
  double mean = gsl_stats_imean (data, n);
  return gsl_stats_iest_stddev_with_mean (data, n, mean) ;
}


double 
gsl_stats_ivariance_with_mean (const int data[], unsigned int n, double mean)
{
  double sum = isum_of_squares (data, n, mean);
  double variance = sum / n;

  return variance;
}

double 
gsl_stats_iest_variance_with_mean (const int data[], unsigned int n,
				  double mean)
{
  double sum = isum_of_squares (data, n, mean);
  double est_variance = sum / (n - 1);

  return est_variance;
}

double 
gsl_stats_istddev_with_mean (const int data[], unsigned int n, double mean)
{
  double sum = isum_of_squares (data, n, mean);
  double stddev = sqrt (sum / n);

  return stddev;
}

double 
gsl_stats_iest_stddev_with_mean (const int data[], unsigned int n,
				double mean)
{
  double sum = isum_of_squares (data, n, mean);
  double est_stddev = sqrt (sum / (n - 1));

  return est_stddev;
}

double 
isum_of_squares (const int data[], unsigned int n, double mean)
{
  /* takes an integer dataset and finds the variance */

  double sum = 0, sum2 = 0;
  unsigned int i;

  /* find the sum of the squares */
  for (i = 0; i < n; i++)
    {
      double delta = (data[i] - mean);
      sum += delta ;  /* FIXME: round off correction, does it work??*/
      sum2 += delta * delta;
    }

  return sum2 - (sum * sum / n) ;
}

