#include <stdlib.h>
#include <gsl_statistics_int.h>

typedef int cmp_fn_t(const void *, const void *) ;

int compare_ints (const int * x, const int * y) ;

int
gsl_stats_int_sort_data (int data[], unsigned int n)
{  
  qsort(data, n, sizeof(int), (cmp_fn_t *) compare_ints) ;
  return 0 ;
}

int
compare_ints (const int * x, const int * y)
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

