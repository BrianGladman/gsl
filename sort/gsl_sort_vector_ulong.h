#ifndef __GSL_SORT_VECTOR_ULONG_H__
#define __GSL_SORT_VECTOR_ULONG_H__

#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_vector_ulong.h>

void gsl_sort_vector_ulong (gsl_vector_ulong * v);
int gsl_sort_vector_ulong_index (gsl_permutation * p, const gsl_vector_ulong * v);

#endif /* __GSL_SORT_VECTOR_ULONG_H__ */
