#ifndef GSL_SORT_VECTOR_CHAR_H
#define GSL_SORT_VECTOR_CHAR_H

#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_vector_char.h>

void gsl_sort_vector_char (gsl_vector_char * v);
void gsl_sort_vector_char_index (gsl_permutation * p, const gsl_vector_char * v);

#endif /* GSL_SORT_VECTOR_CHAR_H */
