
double 
FUNCTION(gsl_stats,lag1_autocorrelation) (const BASE data[], const size_t n)
{
  const double mean = FUNCTION(gsl_stats,mean) (data, n);
  return FUNCTION(gsl_stats,lag1_autocorrelation_with_mean)(data, n, mean);
}

double 
FUNCTION(gsl_stats,lag1_autocorrelation_with_mean) 
     (const BASE data[], const size_t size, const double mean)
{
  /* Compute the lag-1 autocorrelation of a dataset using the
     recurrence relation */

  size_t i;

  long double r1 ;
  long double q = 0 ;
  long double v = (data[0] - mean) * (data[0] - mean) ;

  for (i = 1; i < size ; i++)
    {
      const long double delta0 = (data[i-1] - mean);
      const long double delta1 = (data[i] - mean);
      q += (delta0 * delta1 - q)/(i + 1);
      v += (delta1 * delta1 - v)/(i + 1);
    }

  r1 = q / v ;

  return r1;
}
