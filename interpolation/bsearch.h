/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_INTERP_BSEARCH_H_
#define GSL_INTERP_BSEARCH_H_


unsigned long
interp_bsearch(const double x_array[], double x,
               unsigned long index_lo,
               unsigned long index_hi
               );

#ifdef HAVE_INLINE
extern
inline
unsigned long
interp_bsearch(const double x_array[], double x,
               unsigned long index_lo,
               unsigned long index_hi
               )
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
#endif /* HAVE_INLINE */


#endif /* !GSL_INTERP_BSEARCH_H_ */
