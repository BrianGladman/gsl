#ifndef __GSL_PERMUTE_VECTOR_INT_H__
#define __GSL_PERMUTE_VECTOR_INT_H__

#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_vector_int.h>

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

int gsl_permute_vector_int (const gsl_permutation * p, gsl_vector_int * v);
int gsl_permute_vector_int_inverse (const gsl_permutation * p, gsl_vector_int * v);

__END_DECLS

#endif /* __GSL_PERMUTE_VECTOR_INT_H__ */
