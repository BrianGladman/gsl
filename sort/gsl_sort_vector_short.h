#ifndef __GSL_SORT_VECTOR_SHORT_H__
#define __GSL_SORT_VECTOR_SHORT_H__

#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_vector_short.h>

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

void gsl_sort_vector_short (gsl_vector_short * v);
int gsl_sort_vector_short_index (gsl_permutation * p, const gsl_vector_short * v);

__END_DECLS

#endif /* __GSL_SORT_VECTOR_SHORT_H__ */
