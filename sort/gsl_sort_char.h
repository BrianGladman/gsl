#ifndef __GSL_SORT_CHAR_H__
#define __GSL_SORT_CHAR_H__

#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_permutation.h>

void gsl_sort_char (char * data, size_t stride, size_t n);
int gsl_sort_char_index (size_t * p, const char * data, size_t stride, size_t n);

#endif /* __GSL_SORT_CHAR_H__ */
