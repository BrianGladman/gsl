#ifndef GSL_SORT_VECTOR_UCHAR_H
#define GSL_SORT_VECTOR_UCHAR_H

#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_vector_uchar.h>

void gsl_sort_vector_uchar (gsl_vector_uchar * v);
int gsl_sort_vector_uchar_index (gsl_permutation * p, const gsl_vector_uchar * v);

#endif /* GSL_SORT_VECTOR_UCHAR_H */
