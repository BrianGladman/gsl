#include <math.h>

#include "source.h"

static
double sum_of_squares (const BASE data[], const size_t n,
		       const double mean);

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
  const double sum = sum_of_squares (data, n, mean);
  const double variance = sum / n;

  return variance;
}

double 
FUNCTION(gsl_stats,est_variance_with_mean) (const BASE data[], const size_t n,
					    const double mean)
{
  const double sum = sum_of_squares (data, n, mean);
  const double est_variance = sum / (n - 1);
  
  return est_variance;
}

double 
FUNCTION(gsl_stats,sd_with_mean) (const BASE data[], const size_t n,
				  const double mean)
{
  const double sum = sum_of_squares (data, n, mean);
  const double sd = sqrt (sum / n);

  return sd;
}

double 
FUNCTION(gsl_stats,est_sd_with_mean) (const BASE data[], const size_t n,
				      const double mean)
{
  const double sum = sum_of_squares (data, n, mean);
  const double est_sd = sqrt (sum / (n - 1));

  return est_sd;
}

static double 
sum_of_squares (const BASE data[], const size_t n,
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

