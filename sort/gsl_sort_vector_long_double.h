#ifndef __GSL_SORT_VECTOR_LONG_DOUBLE_H__
#define __GSL_SORT_VECTOR_LONG_DOUBLE_H__

#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_vector_long_double.h>

void gsl_sort_vector_long_double (gsl_vector_long_double * v);
int gsl_sort_vector_long_double_index (gsl_permutation * p, const gsl_vector_long_double * v);

#endif /* __GSL_SORT_VECTOR_LONG_DOUBLE_H__ */
