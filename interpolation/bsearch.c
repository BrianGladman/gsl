#include "bsearch.h"


int
bsearch(const double x_array[], double x, int index_lo, int index_hi)
{
  int ilo = index_lo;
  int ihi = index_hi;
  while(ihi - ilo > 1) {
    int i = (ihi + ilo)/2;
    if(x_array[i] > x)
      ihi = i;
    else
      ilo = i;
  }
  
  return ilo;
}
