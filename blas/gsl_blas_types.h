/*
 * Author:  G. Jungman
 * RCS:     $Id$
 */
/* Based on draft BLAST C interface specification  [Jul 7 1998]
 */
#ifndef GSL_BLAS_TYPES_H_
#define GSL_BLAS_TYPES_H_

#include <sys/types.h>
#include <gsl_complex.h>

typedef  size_t  CBLAS_INDEX;
typedef  size_t  CBLAS_INDEX_t;

enum CBLAS_ORDER      { CblasRowMajor=101, CblasColMajor=102 };
enum CBLAS_TRANSPOSE  { CblasNoTrans=111,  CblasTrans=112,   CblasConjTrans=113 };
enum CBLAS_UPLO       { CblasUpper=121,    CblasLower=122 };
enum CBLAS_DIAG       { CblasNonUnit=131,  CblasUnit=132  };
enum CBLAS_SIDE       { CblasLeft=141,     CblasRight=142 };

typedef  enum CBLAS_ORDER       CBLAS_ORDER_t;
typedef  enum CBLAS_TRANSPOSE   CBLAS_TRANSPOSE_t;
typedef  enum CBLAS_UPLO        CBLAS_UPLO_t;
typedef  enum CBLAS_DIAG        CBLAS_DIAG_t;
typedef  enum CBLAS_SIDE        CBLAS_SIDE_t;

typedef  gsl_complex  COMPLEX;


#endif  /* !GSL_BLAS_TYPES_H_ */
