#ifndef __GSL_SORT_LONG_DOUBLE_H__
#define __GSL_SORT_LONG_DOUBLE_H__

#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_permutation.h>

void gsl_sort_long_double (long double * data, size_t stride, size_t n);
int gsl_sort_long_double_index (size_t * p, const long double * data, size_t stride, size_t n);

#endif /* __GSL_SORT_LONG_DOUBLE_H__ */
