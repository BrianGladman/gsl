/*
 * Author:  G. Jungman
 * RCS:     $Id$
 */
#include "matrix_access.h"
#define ACCESS_UP_CR M_STANDARD_ACCESS_CR
#define ACCESS_UP_CI M_STANDARD_ACCESS_CI
#define ACCESS_LO_CR M_STANDARD_ACCESS_CR
#define ACCESS_LO_CI M_STANDARD_ACCESS_CI
#define KBAND K
#define LDA lda
#include "source_tXsv_c.h"
#undef ACCESS_UP_CR
#undef ACCESS_UP_CI
#undef ACCESS_LO_CR
#undef ACCESS_LO_CI
#undef KBAND
#undef LDA
