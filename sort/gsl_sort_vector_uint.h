#ifndef GSL_SORT_VECTOR_UINT_H
#define GSL_SORT_VECTOR_UINT_H

#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_vector_uint.h>

void gsl_sort_vector_uint (gsl_vector_uint * v);
int gsl_sort_vector_uint_index (gsl_permutation * p, const gsl_vector_uint * v);

#endif /* GSL_SORT_VECTOR_UINT_H */
