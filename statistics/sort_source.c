#include <stdlib.h>

#include "source.h"

typedef int cmp_fn_t(const void *, const void *) ;

static int FUNCTION(my,compare) (const BASE * x, const BASE * y) ;

int
FUNCTION(gsl_stats,sort_data) (BASE data[], const size_t n)
{  
  qsort(data, n, sizeof(BASE), (cmp_fn_t *) FUNCTION(my,compare)) ;
  return 0 ;
}

static int
FUNCTION(my,compare) (const BASE * x, const BASE * y)
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

