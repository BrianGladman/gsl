#ifndef __GSL_SORT_VECTOR_INT_H__
#define __GSL_SORT_VECTOR_INT_H__

#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_vector_int.h>

void gsl_sort_vector_int (gsl_vector_int * v);
int gsl_sort_vector_int_index (gsl_permutation * p, const gsl_vector_int * v);

#endif /* __GSL_SORT_VECTOR_INT_H__ */
