#include <stdlib.h>

static int FUNCTION(my,compare) (const void * x, const void * y) ;

void
FUNCTION(gsl_stats,sort_data) (BASE data[], const size_t n)
{  
  qsort(data, n, sizeof(BASE), FUNCTION(my,compare)) ;
}

static int
FUNCTION(my,compare) (const void * x, const void * y)
{
  if (* (const BASE *)x < * (const BASE *)y) 
    {
      return -1 ;
    }
  else if (* (const BASE *)x > * (const BASE *)y)
    {
      return +1 ;
    } 
  else 
    {
      return 0 ;
    }
}

