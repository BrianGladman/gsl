#ifndef __GSL_SORT_ULONG_H__
#define __GSL_SORT_ULONG_H__

#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_permutation.h>

void gsl_sort_ulong (unsigned long * data, size_t stride, size_t n);
int gsl_sort_ulong_index (size_t * p, const unsigned long * data, size_t stride, size_t n);

#endif /* __GSL_SORT_ULONG_H__ */
