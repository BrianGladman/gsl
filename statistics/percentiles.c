#include <gsl_statistics.h>

double
gsl_stats_percentile_from_sorted_data (const double sorted_data[], 
				       unsigned int n, double f)
{
  double index = f * (n - 1) ;
  unsigned int lhs = (int)index ;
  double delta = index - lhs ;
  double result;

  if (n == 0)
    return 0.0 ;

  if (lhs == n - 1)
    {
      result = sorted_data[lhs] ;
    }
  else 
    {
      result = (1 - delta) * sorted_data[lhs] + delta * sorted_data[lhs + 1] ;
    }

  return result ;
}
