#include <config.h> 
#include <stdlib.h>
#include "utils.h"

#undef BASE
#include "init_source.c"
#undef BASE

#define BASE int
#include "init_source.c"
#undef BASE
