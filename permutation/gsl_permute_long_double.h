#ifndef __GSL_PERMUTE_LONG_DOUBLE_H__
#define __GSL_PERMUTE_LONG_DOUBLE_H__

#include <stdlib.h>
#include <gsl/gsl_errno.h>
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

int gsl_permute_long_double (const size_t * p, long double * data, size_t stride, size_t n);
int gsl_permute_long_double_inverse (const size_t * p, long double * data, size_t stride, size_t n);

__END_DECLS

#endif /* __GSL_PERMUTE_LONG_DOUBLE_H__ */
