#ifndef __GSL_SORT_VECTOR_LONG_H__
#define __GSL_SORT_VECTOR_LONG_H__

#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_vector_long.h>

void gsl_sort_vector_long (gsl_vector_long * v);
int gsl_sort_vector_long_index (gsl_permutation * p, const gsl_vector_long * v);

#endif /* __GSL_SORT_VECTOR_LONG_H__ */
