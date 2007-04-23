#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_errno.h>

#include "../matrix/view.h"

#define BASE_DOUBLE
#include "templates_on.h"
#include "subrowcol_source.c"
#include "templates_off.h"
#undef  BASE_DOUBLE

#define USE_QUALIFIER
#define QUALIFIER const

#define BASE_DOUBLE
#include "templates_on.h"
#include "subrowcol_source.c"
#include "templates_off.h"
#undef  BASE_DOUBLE
