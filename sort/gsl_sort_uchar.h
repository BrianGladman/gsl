#ifndef GSL_SORT_UCHAR_H
#define GSL_SORT_UCHAR_H

#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_permutation.h>

void gsl_sort_uchar (unsigned char * data, size_t stride, size_t n);
int gsl_sort_uchar_index (size_t * p, const unsigned char * data, size_t stride, size_t n);

#endif /* GSL_SORT_UCHAR_H */
