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
#include "complex_internal.h"
#include "gsl_blas_raw_cm.h"

