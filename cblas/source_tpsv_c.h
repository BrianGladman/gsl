/*
 * Author:  G. Jungman
 * RCS:     $Id$
 */
#include "matrix_access.h"
#define ACCESS_UP_CR M_PACKEDTRUP_ACCESS_CR
#define ACCESS_UP_CI M_PACKEDTRUP_ACCESS_CI
#define ACCESS_LO_CR M_PACKEDTRLO_ACCESS_CR
#define ACCESS_LO_CI M_PACKEDTRLO_ACCESS_CI
#define KBAND (N-1)
#define LDA 0
#include "source_tXsv_c.h"
#undef ACCESS_UP_CR
#undef ACCESS_UP_CI
#undef ACCESS_LO_CR
#undef ACCESS_LO_CI
#undef KBAND
#undef LDA
