/*
 * Author:  G. Jungman
 * RCS:     $Id$
 */
#include "matrix_access.h"
#define ACCESS_UP M_STANDARD_ACCESS
#define ACCESS_LO M_STANDARD_ACCESS
#define KBAND K
#define LDA lda
#include "source_tXsv_r.h"
#undef ACCESS_UP
#undef ACCESS_LO
#undef KBAND
#undef LDA
