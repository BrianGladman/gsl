#ifndef __GSL_HEAPSORT_H__
#define __GSL_HEAPSORT_H__

#include <gsl/gsl_permutation.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

typedef int (*gsl_comparison_fn_t) (const void *, const void *);

void gsl_heapsort (void * array, size_t count, size_t size, gsl_comparison_fn_t compare);
int gsl_heapsort_index (size_t * p, const void * array, size_t count, size_t size, gsl_comparison_fn_t compare);

__END_DECLS

#endif /* __GSL_HEAPSORT_H__ */
