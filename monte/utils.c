#include <config.h> 
#include <stdlib.h>
#include <gsl_errno.h>

#include "utils.h"
#include "source.h"

#undef BASE
#include "init_source.c"
#undef BASE

#define BASE int
#include "init_source.c"
#undef BASE
