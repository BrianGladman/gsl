#ifndef GSL_SORT_INT_H
#define GSL_SORT_INT_H

#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_permutation.h>

void gsl_sort_int (int * data, size_t stride, size_t n);
int gsl_sort_int_index (size_t * p, const int * data, size_t stride, size_t n);

#endif /* GSL_SORT_INT_H */
