#ifndef __GSL_SORT_VECTOR_USHORT_H__
#define __GSL_SORT_VECTOR_USHORT_H__

#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_vector_ushort.h>

void gsl_sort_vector_ushort (gsl_vector_ushort * v);
int gsl_sort_vector_ushort_index (gsl_permutation * p, const gsl_vector_ushort * v);

#endif /* __GSL_SORT_VECTOR_USHORT_H__ */
