#ifndef __GSL_SORT_USHORT_H__
#define __GSL_SORT_USHORT_H__

#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_permutation.h>

void gsl_sort_ushort (unsigned short * data, size_t stride, size_t n);
int gsl_sort_ushort_index (size_t * p, const unsigned short * data, size_t stride, size_t n);

#endif /* __GSL_SORT_USHORT_H__ */
