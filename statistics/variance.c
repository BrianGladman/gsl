#include <math.h>
#include <gsl_statistics.h>

double sum_of_squares (const double data[], const size_t n,
		       const double mean);

double 
gsl_stats_variance (const double data[], const size_t n)
{
  const double mean = gsl_stats_mean (data, n);
  return gsl_stats_variance_with_mean(data, n, mean);
}

double 
gsl_stats_est_variance (const double data[], const size_t n)
{
  const double mean = gsl_stats_mean (data, n);
  return gsl_stats_est_variance_with_mean(data, n, mean);

}

double 
gsl_stats_sd (const double data[], const size_t n)
{
  const double mean = gsl_stats_mean (data, n);
  return gsl_stats_sd_with_mean (data, n, mean) ;
}

double 
gsl_stats_est_sd (const double data[], const size_t n)
{
  const double mean = gsl_stats_mean (data, n);
  return gsl_stats_est_sd_with_mean (data, n, mean) ;
}


double 
gsl_stats_variance_with_mean (const double data[], const size_t n, 
			      const double mean)
{
  const double sum = sum_of_squares (data, n, mean);
  const double variance = sum / n;

  return variance;
}

double 
gsl_stats_est_variance_with_mean (const double data[], const size_t n,
				  const double mean)
{
  const double sum = sum_of_squares (data, n, mean);
  const double est_variance = sum / (n - 1);

  return est_variance;
}

double 
gsl_stats_sd_with_mean (const double data[], const size_t n,
			const double mean)
{
  const double sum = sum_of_squares (data, n, mean);
  const double sd = sqrt (sum / n);

  return sd;
}

double 
gsl_stats_est_sd_with_mean (const double data[], const size_t n,
			    const double mean)
{
  const double sum = sum_of_squares (data, n, mean);
  const double est_sd = sqrt (sum / (n - 1));

  return est_sd;
}

double 
sum_of_squares (const double data[], const size_t n,
		const double mean)
{
  /* takes a dataset and finds the variance */

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

