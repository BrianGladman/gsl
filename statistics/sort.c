#include <stdlib.h>
#include <gsl_statistics.h>

typedef int cmp_fn_t(const void *, const void *) ;

int compare_doubles (const double * x, const double * y) ;

int
gsl_stats_sort_data (double data[], unsigned int n)
{  
  qsort(data, n, sizeof(double), (cmp_fn_t *) compare_doubles) ;
  return 0 ;
}

int
compare_doubles (const double * x, const double * y)
{
  if (*x < *y) 
    {
      return -1 ;
    }
  else if (*x > *y)
    {
      return +1 ;
    } 
  else 
    {
      return 0 ;
    }
}

