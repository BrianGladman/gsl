#ifndef __GSL_SORT_LONG_H__
#define __GSL_SORT_LONG_H__

#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_permutation.h>

void gsl_sort_long (long * data, size_t stride, size_t n);
int gsl_sort_long_index (size_t * p, const long * data, size_t stride, size_t n);

#endif /* __GSL_SORT_LONG_H__ */
