#ifndef __GSL_SORT_SHORT_H__
#define __GSL_SORT_SHORT_H__

#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_permutation.h>

void gsl_sort_short (short * data, size_t stride, size_t n);
int gsl_sort_short_index (size_t * p, const short * data, size_t stride, size_t n);

#endif /* __GSL_SORT_SHORT_H__ */
