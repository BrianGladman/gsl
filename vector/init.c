#include <config.h>
#include <gsl_vector.h>

#define BASE_GSL_COMPLEX_LONG
#define IN_FORMAT "%lg %lg"
#include "source.h"
#include "init_source.c"
#include "unsource.h"
#undef  IN_FORMAT
#undef  BASE_GSL_COMPLEX_LONG

#define BASE_GSL_COMPLEX
#define IN_FORMAT "%lg %lg"
#include "source.h"
#include "init_source.c"
#include "unsource.h"
#undef  IN_FORMAT
#undef  BASE_GSL_COMPLEX

#define BASE_GSL_COMPLEX_FLOAT
#define IN_FORMAT "%g %g"
#include "source.h"
#include "init_source.c"
#include "unsource.h"
#undef  IN_FORMAT
#undef  BASE_GSL_COMPLEX_FLOAT

#define BASE_LONG_DOUBLE
#define IN_FORMAT "%lg"
#include "source.h"
#include "init_source.c"
#include "unsource.h"
#undef  IN_FORMAT
#undef  BASE_LONG_DOUBLE

#define BASE_DOUBLE
#define IN_FORMAT "%lg"
#include "source.h"
#include "init_source.c"
#include "unsource.h"
#undef  IN_FORMAT
#undef  BASE_DOUBLE

#define BASE_FLOAT
#define IN_FORMAT "%lg"
#include "source.h"
#include "init_source.c"
#include "unsource.h"
#undef  IN_FORMAT
#undef  BASE_FLOAT

#define BASE_ULONG
#define IN_FORMAT "%lg"
#include "source.h"
#include "init_source.c"
#include "unsource.h"
#undef  IN_FORMAT
#undef  BASE_ULONG

#define BASE_LONG
#define IN_FORMAT "%lg"
#include "source.h"
#include "init_source.c"
#include "unsource.h"
#undef  IN_FORMAT
#undef  BASE_LONG

#define BASE_UINT
#define IN_FORMAT "%ld"
#include "source.h"
#include "init_source.c"
#include "unsource.h"
#undef  IN_FORMAT
#undef  BASE_UINT

#define BASE_INT
#define IN_FORMAT "%ld"
#include "source.h"
#include "init_source.c"
#include "unsource.h"
#undef  IN_FORMAT
#undef  BASE_INT

#define BASE_USHORT
#define IN_FORMAT "%ld"
#include "source.h"
#include "init_source.c"
#include "unsource.h"
#undef  IN_FORMAT
#undef  BASE_USHORT

#define BASE_SHORT
#define IN_FORMAT "%ld"
#include "source.h"
#include "init_source.c"
#include "unsource.h"
#undef  IN_FORMAT
#undef  BASE_SHORT

#define BASE_UCHAR
#define IN_FORMAT "%ld"
#include "source.h"
#include "init_source.c"
#include "unsource.h"
#undef  IN_FORMAT
#undef  BASE_UCHAR

#define BASE_CHAR
#define IN_FORMAT "%ld"
#include "source.h"
#include "init_source.c"
#include "unsource.h"
#undef  IN_FORMAT
#undef  BASE_CHAR
