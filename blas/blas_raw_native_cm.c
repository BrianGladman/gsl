/*
 * Author:  B. Gough
 * RCS:     $Id$
 */
/* Native implementation of column major operations.
 * Conforms to gsl_blas_raw_cm interface.
 * This includes only the level 2 and level 3 functions,
 * since the level 1 functions are insensitive to the
 * storage scheme.
 */
#include <math.h>
#include "gsl_blas_raw_cm.h"



/* ===========================================================================
 * Generate the col-major versions of the level 2 and level 3 operations.
 * ===========================================================================
 */

#define  MACCESS(s, i, j)  ((s)*(j) + (i))
#define  FUNC(f) f ## _ ## cm
#include "blas_raw_native_L23_source.c"
#undef   FUNC
#undef   MACCESS

