#ifndef GSL_SORT_H
#define GSL_SORT_H

typedef int (*gsl_comparison_fn_t) (const void *, const void *);

void gsl_sort (void * array, size_t count, size_t size, gsl_comparison_fn_t compare);

#endif /* GSL_SORT_H */
