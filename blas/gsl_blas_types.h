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
typedef  enum { CblasRowMajor=101, CblasColMajor=102 }                    CBLAS_ORDER;
typedef  enum { CblasNoTrans=111,  CblasTrans=112,   CblasConjTrans=113 } CBLAS_TRANSPOSE;
typedef  enum { CblasUpper=121,    CblasLower=122 }                       CBLAS_UPLO;
typedef  enum { CblasNonUnit=131,  CblasUnit=132 }                        CBLAS_DIAG;
typedef  enum { CblasLeft=141,     CblasRight=142 }                       CBLAS_SIDE;

typedef  gsl_complex  COMPLEX;


#endif  /* !GSL_BLAS_TYPES_H_ */
