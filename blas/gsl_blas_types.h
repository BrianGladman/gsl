/*
 * Author:  G. Jungman
 * RCS:     $Id$
 */
/* Based on draft BLAST C interface specification  [Jul 7 1998]
 */
#ifndef __GSL_BLAS_TYPES_H__
#define __GSL_BLAS_TYPES_H__

#include <sys/types.h>
#include <gsl/gsl_complex.h>

#include <gsl/gsl_cblas.h>

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

typedef  CBLAS_INDEX  CBLAS_INDEX_t;
typedef  enum CBLAS_ORDER       CBLAS_ORDER_t;
typedef  enum CBLAS_TRANSPOSE   CBLAS_TRANSPOSE_t;
typedef  enum CBLAS_UPLO        CBLAS_UPLO_t;
typedef  enum CBLAS_DIAG        CBLAS_DIAG_t;
typedef  enum CBLAS_SIDE        CBLAS_SIDE_t;

/* typedef  gsl_complex  COMPLEX; */

__END_DECLS


#endif /* __GSL_BLAS_TYPES_H__ */
