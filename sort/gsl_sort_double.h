#ifndef __GSL_SORT_DOUBLE_H__
#define __GSL_SORT_DOUBLE_H__

#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_permutation.h>

void gsl_sort (double * data, size_t stride, size_t n);
int gsl_sort_index (size_t * p, const double * data, size_t stride, size_t n);

#endif /* __GSL_SORT_DOUBLE_H__ */
