#ifndef __GSL_PERMUTE_VECTOR_FLOAT_H__
#define __GSL_PERMUTE_VECTOR_FLOAT_H__

#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_vector_float.h>

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

int gsl_permute_vector_float (const gsl_permutation * p, gsl_vector_float * v);
int gsl_permute_vector_float_inverse (const gsl_permutation * p, gsl_vector_float * v);

__END_DECLS

#endif /* __GSL_PERMUTE_VECTOR_FLOAT_H__ */
