#include <math.h>
#include <gsl_statistics.h>

double sum_of_squares (const double data[], unsigned int n, double mean);

double 
gsl_stats_variance (const double data[], unsigned int n)
{
  double mean = gsl_stats_mean (data, n);
  return gsl_stats_variance_with_mean(data, n, mean);
}

double 
gsl_stats_est_variance (const double data[], unsigned int n)
{
  double mean = gsl_stats_mean (data, n);
  return gsl_stats_est_variance_with_mean(data, n, mean);

}

double 
gsl_stats_sd (const double data[], unsigned int n)
{
  double mean = gsl_stats_mean (data, n);
  return gsl_stats_sd_with_mean (data, n, mean) ;
}

double 
gsl_stats_est_sd (const double data[], unsigned int n)
{
  double mean = gsl_stats_mean (data, n);
  return gsl_stats_est_sd_with_mean (data, n, mean) ;
}


double 
gsl_stats_variance_with_mean (const double data[], unsigned int n, double mean)
{
  double sum = sum_of_squares (data, n, mean);
  double variance = sum / n;

  return variance;
}

double 
gsl_stats_est_variance_with_mean (const double data[], unsigned int n,
				  double mean)
{
  double sum = sum_of_squares (data, n, mean);
  double est_variance = sum / (n - 1);

  return est_variance;
}

double 
gsl_stats_sd_with_mean (const double data[], unsigned int n, double mean)
{
  double sum = sum_of_squares (data, n, mean);
  double sd = sqrt (sum / n);

  return sd;
}

double 
gsl_stats_est_sd_with_mean (const double data[], unsigned int n,
				double mean)
{
  double sum = sum_of_squares (data, n, mean);
  double est_sd = sqrt (sum / (n - 1));

  return est_sd;
}

double 
sum_of_squares (const double data[], unsigned int n, double mean)
{
  /* takes a dataset and finds the variance */

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

