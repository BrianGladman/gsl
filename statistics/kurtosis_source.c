double 
FUNCTION(gsl_stats,kurtosis) (const BASE data[], const size_t n)
{
  const double mean = FUNCTION(gsl_stats,mean)(data, n);
  const double est_sd = FUNCTION(gsl_stats,est_sd_with_mean)(data, n, mean);
  return FUNCTION(gsl_stats,kurtosis_with_mean_and_sd)(data, n, mean, est_sd);
}
    
double 
FUNCTION(gsl_stats,kurtosis_with_mean_and_sd) (const BASE data[], 
					       const size_t n,
					       const double mean, 
					       const double sd)
{
  /* takes a dataset and finds the kurtosis */

  long double avg = 0, kurtosis;
  size_t i;

  /* find the fourth moment the deviations, normalized by the sd */

  /* we use a recurrence relation to stably update a running value so
     there aren't any large sums that can overflow */

  for (i = 0; i < n; i++)
    {
      const long double x = (data[i] - mean) / sd;
      avg += (x * x * x * x - avg)/(i + 1);
    }

  kurtosis = avg - 3.0;  /* makes kurtosis zero for a Gaussian */

  return kurtosis;
}

