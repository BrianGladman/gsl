
double
FUNCTION(gsl_stats,quantile_from_sorted_data) (const BASE sorted_data[], 
						 const size_t n,
						 const double f)
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
