/*
 * Author:  G. Jungman
 * RCS:     $Id$
 */
/* Raw BLAS interface for row-major matrices.
 * Based on draft BLAST C interface specification  [Jul 7 1998]
 */
#ifndef GSL_BLAS_RAW_H
#define GSL_BLAS_RAW_H

#include <gsl/gsl_blas_raw_L1.h>
#include <gsl/gsl_blas_raw_L2.h>
#include <gsl/gsl_blas_raw_L3.h>


#if defined(HAVE_INLINE) && defined(HAVE_CBLAS)
#include <cblas.h>

/* insert inline cblas implementation of above here */

#endif /* defined(HAVE_INLINE) && defined(HAVE_CBLAS) */


#endif /* !GSL_BLAS_RAW_H */
