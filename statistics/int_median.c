/* This is an automatically generated file created by
   convert_double_to_int.pl -- do not edit                */

#include <gsl_statistics_int.h>

double
gsl_stats_int_median_from_sorted_data (const int sorted_data[], const size_t n)
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

