#ifndef __GSL_SORT_VECTOR_FLOAT_H__
#define __GSL_SORT_VECTOR_FLOAT_H__

#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_vector_float.h>

void gsl_sort_vector_float (gsl_vector_float * v);
int gsl_sort_vector_float_index (gsl_permutation * p, const gsl_vector_float * v);

#endif /* __GSL_SORT_VECTOR_FLOAT_H__ */
