#ifndef __GSL_SORT_INT_H__
#define __GSL_SORT_INT_H__

#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_permutation.h>

void gsl_sort_int (int * data, size_t stride, size_t n);
int gsl_sort_int_index (size_t * p, const int * data, size_t stride, size_t n);

#endif /* __GSL_SORT_INT_H__ */
