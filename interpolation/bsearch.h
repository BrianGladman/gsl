#ifndef GSL_BSEARCH_H_
#define GSL_BSEARCH_H_

int
bsearch(const double x_array[], double x, int index_lo, int index_hi);

#define CHECK_BSEARCH(xa, x, ilo, ihi) ( (ihi)-(ilo) == 1 ? (ilo) : bsearch(xa, x, ilo, ihi) )

#endif  /* !GSL_BSEARCH_H_ */
