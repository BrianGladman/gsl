#include <config.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <stdio.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_complex.h>

#include <gsl/gsl_fft_complex.h>

#define BASE_DOUBLE
#include "templates_on.h"
#include "compare_source.c"
#include "templates_off.h"
#undef  BASE_DOUBLE

#define BASE_FLOAT
#include "templates_on.h"
#include "compare_source.c"
#include "templates_off.h"
#undef  BASE_FLOAT
