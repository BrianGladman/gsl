#ifndef __GSL_PERMUTE_UINT_H__
#define __GSL_PERMUTE_UINT_H__

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

int gsl_permute_uint (const size_t * p, unsigned int * data, size_t stride, size_t n);
int gsl_permute_uint_inverse (const size_t * p, unsigned int * data, size_t stride, size_t n);

__END_DECLS

#endif /* __GSL_PERMUTE_UINT_H__ */
