double
FUNCTION(gsl_stats,absdev) (const BASE data[], const size_t n)
{
  const double mean = FUNCTION(gsl_stats,mean)(data, n);
  return FUNCTION(gsl_stats,absdev_with_mean)(data, n, mean);
}
    
double 
FUNCTION(gsl_stats,absdev_with_mean) (const BASE data[], 
				      const size_t n, 
				      const double mean)
{
  /* takes a dataset and finds the absolute deviation */

  double sum = 0, absdev;
  size_t i;

  /* find the sum of the absolute deviations */
  for (i = 0; i < n; i++)
    {
      const double delta = fabs(data[i] - mean);
      sum += delta;
    }

  absdev = sum / n;

  return absdev;
}

