/* This is an automatically generated file created by
   convert_double_to_int.pl -- do not edit                */

#include <math.h>
#include <gsl_statistics_int.h>

double int_sum_of_squares (const int data[], const size_t n,
		       const double mean);

double 
gsl_stats_int_variance (const int data[], const size_t n)
{
  const double mean = gsl_stats_int_mean (data, n);
  return gsl_stats_int_variance_with_mean(data, n, mean);
}

double 
gsl_stats_int_est_variance (const int data[], const size_t n)
{
  const double mean = gsl_stats_int_mean (data, n);
  return gsl_stats_int_est_variance_with_mean(data, n, mean);

}

double 
gsl_stats_int_sd (const int data[], const size_t n)
{
  const double mean = gsl_stats_int_mean (data, n);
  return gsl_stats_int_sd_with_mean (data, n, mean) ;
}

double 
gsl_stats_int_est_sd (const int data[], const size_t n)
{
  const double mean = gsl_stats_int_mean (data, n);
  return gsl_stats_int_est_sd_with_mean (data, n, mean) ;
}


double 
gsl_stats_int_variance_with_mean (const int data[], const size_t n, 
			      const double mean)
{
  const double sum = int_sum_of_squares (data, n, mean);
  const double variance = sum / n;

  return variance;
}

double 
gsl_stats_int_est_variance_with_mean (const int data[], const size_t n,
				  const double mean)
{
  const double sum = int_sum_of_squares (data, n, mean);
  const double est_variance = sum / (n - 1);

  return est_variance;
}

double 
gsl_stats_int_sd_with_mean (const int data[], const size_t n,
			const double mean)
{
  const double sum = int_sum_of_squares (data, n, mean);
  const double sd = sqrt (sum / n);

  return sd;
}

double 
gsl_stats_int_est_sd_with_mean (const int data[], const size_t n,
			    const double mean)
{
  const double sum = int_sum_of_squares (data, n, mean);
  const double est_sd = sqrt (sum / (n - 1));

  return est_sd;
}

double 
int_sum_of_squares (const int data[], const size_t n,
		const double mean)
{
  /* takes an integer dataset and finds the variance */

  double sum = 0, sum2 = 0;
  size_t i;

  /* find the sum of the squares */
  for (i = 0; i < n; i++)
    {
      const double delta = (data[i] - mean);
      sum += delta ;  /* FIXME: round off correction, does it work??*/
      sum2 += delta * delta;
    }

  return sum2 - (sum * sum / n) ;
}

