#include <math.h>

#include "source.h"

static double 
compute_variance (const BASE data[], const size_t n, const double mean);

double 
FUNCTION(gsl_stats,variance) (const BASE data[], const size_t n)
{
  const double mean = FUNCTION(gsl_stats,mean) (data, n);
  return FUNCTION(gsl_stats,variance_with_mean)(data, n, mean);
}

double 
FUNCTION(gsl_stats,est_variance) (const BASE data[], const size_t n)
{
  const double mean = FUNCTION(gsl_stats,mean) (data, n);
  return FUNCTION(gsl_stats,est_variance_with_mean)(data, n, mean);

}

double 
FUNCTION(gsl_stats,sd) (const BASE data[], const size_t n)
{
  const double mean = FUNCTION(gsl_stats,mean) (data, n);
  return FUNCTION(gsl_stats,sd_with_mean) (data, n, mean) ;
}

double 
FUNCTION(gsl_stats,est_sd) (const BASE data[], const size_t n)
{
  const double mean = FUNCTION(gsl_stats,mean) (data, n);
  return FUNCTION(gsl_stats,est_sd_with_mean) (data, n, mean) ;
}


double 
FUNCTION(gsl_stats,variance_with_mean) (const BASE data[], const size_t n, 
					const double mean)
{
  const double variance = compute_variance (data, n, mean);
  return variance;
}

double 
FUNCTION(gsl_stats,est_variance_with_mean) (const BASE data[], const size_t n,
					    const double mean)
{
  const double variance = compute_variance (data, n, mean);
  const double est_variance = variance * ((double)n / (double)(n - 1));
  
  return est_variance;
}

double 
FUNCTION(gsl_stats,sd_with_mean) (const BASE data[], const size_t n,
				  const double mean)
{
  const double variance = compute_variance (data, n, mean);
  const double sd = sqrt (variance);

  return sd;
}

double 
FUNCTION(gsl_stats,est_sd_with_mean) (const BASE data[], const size_t n,
				      const double mean)
{
  const double variance = compute_variance (data, n, mean);
  const double est_sd = sqrt (variance * ((double)n / (double)(n - 1)));

  return est_sd;
}

static double 
compute_variance (const BASE data[], const size_t n, const double mean)
{
  /* takes a dataset and finds the variance */

  double variance = 0 ;

  size_t i;

  /* find the sum of the squares */
  for (i = 0; i < n; i++)
    {
      const double delta = (data[i] - mean);
      variance += (delta * delta - variance) / ((double)(i + 1));
    }

  return variance ;
}

