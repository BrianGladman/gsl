#include <gsl_statistics.h>

double
gsl_stats_median_from_sorted_data (const double sorted_data[], const size_t n)
{
  double median ;
  const size_t lhs = (n - 1) / 2 ;
  const size_t rhs = n / 2 ;
  
  if (n == 0)
    return 0.0 ;

  if (lhs == rhs)
    {
      median = sorted_data[lhs] ;
    }
  else 
    {
      median = (sorted_data[lhs] + sorted_data[rhs])/2.0 ;
    }

  return median ;
}

