#ifndef __GSL_SORT_VECTOR_DOUBLE_H__
#define __GSL_SORT_VECTOR_DOUBLE_H__

#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_vector_double.h>

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

void gsl_sort_vector (gsl_vector * v);
int gsl_sort_vector_index (gsl_permutation * p, const gsl_vector * v);

__END_DECLS

#endif /* __GSL_SORT_VECTOR_DOUBLE_H__ */
