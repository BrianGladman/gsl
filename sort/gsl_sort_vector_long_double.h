#ifndef GSL_SORT_VECTOR_LONG_DOUBLE_H
#define GSL_SORT_VECTOR_LONG_DOUBLE_H

#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_vector_long_double.h>

void gsl_sort_vector_long_double (gsl_vector_long_double * v);
void gsl_sort_vector_long_double_index (gsl_permutation * p, const gsl_vector_long_double * v);

#endif /* GSL_SORT_VECTOR_LONG_DOUBLE_H */
