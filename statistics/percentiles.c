#include <gsl_statistics.h>

double
gsl_stats_percentile_from_sorted_data (const double sorted_data[], 
				       const size_t n, const double f)
{
  const double index = f * (n - 1) ;
  const size_t lhs = (int)index ;
  const double delta = index - lhs ;
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
