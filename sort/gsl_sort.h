#ifndef GSL_SORT_H
#define GSL_SORT_H

#include <gsl/gsl_permutation.h>

typedef int (*gsl_comparison_fn_t) (const void *, const void *);

void gsl_sort (void * array, size_t count, size_t size, gsl_comparison_fn_t compare);
int gsl_sort_index (gsl_permutation * permutation, const void * array, size_t count, size_t size, gsl_comparison_fn_t compare);

#endif /* GSL_SORT_H */
