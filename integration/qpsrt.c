#include <config.h>
#include <stdlib.h>

#include "qpsrt.h"

void qpsrt (const size_t limit, const size_t last, 
	    size_t * maxerr_index, double * maxerr_value, 
	    const double elist[], size_t iord[], 
	    size_t * nrmax)
{
  double errmax ;
  double errmin ;
  int i, k, top;

  size_t i_nrmax = * nrmax ;
  size_t i_maxerr = iord[i_nrmax] ;
  
  /* Check whether the list contains more than two error estimates */

  if (last < 2) 
    {
      iord[0] = 0 ;
      iord[1] = 1 ;
      * maxerr_index = i_maxerr ;
      * maxerr_value = elist[i_maxerr] ;
      return ;
    }

  errmax = elist[i_maxerr] ;

  /* This part of the routine is only executed if, due to a difficult
     integrand, subdivision increased the error estimate. In the normal
     case the insert procedure should start after the nrmax-th largest
     error estimate. */

  while (i_nrmax > 0 && errmax > elist[iord[i_nrmax - 1]]) 
    {
      iord[i_nrmax] = iord[i_nrmax - 1] ;
      /* printf("NR copied iord[%d] to iord[%d]\n", i_nrmax-1, i_nrmax); */
      i_nrmax-- ;
    } 

  /* Compute the number of elements in the list to be maintained in
     descending order. This number depends on the number of
     subdivisions still allowed. */
  
  if(last < (limit/2 + 2)) 
    {
      top = last ;
    }
  else
    {
      top = limit - last + 1;
    }

  /* printf("top = %d\n", top) ; */

  
  /* Insert errmax by traversing the list top-down, starting
     comparison from the element elist(iord(i_nrmax+1)). */
  
  i = i_nrmax + 1 ;
  
  while (errmax < elist[iord[i]] && i < top)
    {
      iord[i-1] = iord[i] ;
      /* printf("DOWN copied iord[%d] to iord[%d]\n", i, i-1); */
      i++ ;
    }
  
  iord[i-1] = i_maxerr ;
  /* printf("I-1 put %d into iord[%d]\n", i_maxerr, i-1) ; */
  
  /* Insert errmin by traversing the list bottom-up */
  
  errmin = elist[last] ;
  
  k = top - 1 ;
  
  while (errmin >= elist[iord[k]] && k > i - 2)
    {
      iord[k+1] = iord[k] ;
      /* printf("UP copied iord[%d] to iord[%d]\n", k, k+1); */
      k-- ;
    }
  
  iord[k+1] = last ;
  /* printf("K put %d into iord[%d]\n", last, k+1) ; */

  /* Set maxerr_index and maxerr_value */

  i_maxerr = iord[i_nrmax] ;
  
  * maxerr_index = i_maxerr ;
  * maxerr_value = elist[i_maxerr] ;

  * nrmax = i_nrmax ;
}

