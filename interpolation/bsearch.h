/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_INTERP_BSEARCH_H_
#define GSL_INTERP_BSEARCH_H_

size_t
gsl_interp_bsearch(const double x_array[], double x,
               size_t index_lo,
               size_t index_hi
               );

#ifdef HAVE_INLINE
extern
inline
size_t
gsl_interp_bsearch(const double x_array[], double x,
               size_t index_lo,
               size_t index_hi
               )
{
  size_t ilo = index_lo;
  size_t ihi = index_hi;
  while(ihi > ilo + 1) {
    size_t i = (ihi + ilo)/2;
    if(x_array[i] > x)
      ihi = i;
    else
      ilo = i;
  }
  
  return ilo;
}
#endif /* HAVE_INLINE */


#endif /* !GSL_INTERP_BSEARCH_H_ */
