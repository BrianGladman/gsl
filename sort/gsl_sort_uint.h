#ifndef __GSL_SORT_UINT_H__
#define __GSL_SORT_UINT_H__

#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_permutation.h>

void gsl_sort_uint (unsigned int * data, size_t stride, size_t n);
int gsl_sort_uint_index (size_t * p, const unsigned int * data, size_t stride, size_t n);

#endif /* __GSL_SORT_UINT_H__ */
