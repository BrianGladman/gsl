#ifndef __GSL_SORT_VECTOR_SHORT_H__
#define __GSL_SORT_VECTOR_SHORT_H__

#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_vector_short.h>

void gsl_sort_vector_short (gsl_vector_short * v);
int gsl_sort_vector_short_index (gsl_permutation * p, const gsl_vector_short * v);

#endif /* __GSL_SORT_VECTOR_SHORT_H__ */
