/*
 * Author:  G. Jungman
 * RCS:     $Id$
 */
#include "matrix_access.h"
#define ACCESS_UP M_PACKEDTRUP_ACCESS
#define ACCESS_LO M_PACKEDTRLO_ACCESS
#define KBAND (N-1)
#include "source_tXsv_r.h"
#undef ACCESS_UP
#undef ACCESS_LO
#undef KBAND
