#include <config.h> 
#include <stdlib.h>
#include <gsl_errno.h>

#include "utils.h"

#define BASE_DOUBLE
#include "templates_on.h"
#include "init_source.c"
#include "templates_off.h"
#undef BASE_DOUBLE

#define BASE_INT
#include "templates_on.h"
#include "init_source.c"
#include "templates_off.h"
#undef BASE_INT
