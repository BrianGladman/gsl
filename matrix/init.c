#include <config.h>
#include <gsl_matrix.h>

#define BASE_GSL_COMPLEX_LONG
#include "source.h"
#include "init_source.c"
#include "unsource.h"
#undef  BASE_GSL_COMPLEX_LONG

#define BASE_GSL_COMPLEX
#include "source.h"
#include "init_source.c"
#include "unsource.h"
#undef  BASE_GSL_COMPLEX

#define BASE_GSL_COMPLEX_FLOAT
#include "source.h"
#include "init_source.c"
#include "unsource.h"
#undef  BASE_GSL_COMPLEX_FLOAT

#define BASE_LONG_DOUBLE
#include "source.h"
#include "init_source.c"
#include "unsource.h"
#undef  BASE_LONG_DOUBLE

#define BASE_DOUBLE
#include "source.h"
#include "init_source.c"
#include "unsource.h"
#undef  BASE_DOUBLE

#define BASE_FLOAT
#include "source.h"
#include "init_source.c"
#include "unsource.h"
#undef  BASE_FLOAT

#define BASE_ULONG
#include "source.h"
#include "init_source.c"
#include "unsource.h"
#undef  BASE_ULONG

#define BASE_LONG
#include "source.h"
#include "init_source.c"
#include "unsource.h"
#undef  BASE_LONG

#define BASE_UINT
#include "source.h"
#include "init_source.c"
#include "unsource.h"
#undef  BASE_UINT

#define BASE_INT
#include "source.h"
#include "init_source.c"
#include "unsource.h"
#undef  BASE_INT

#define BASE_USHORT
#include "source.h"
#include "init_source.c"
#include "unsource.h"
#undef  BASE_USHORT

#define BASE_SHORT
#include "source.h"
#include "init_source.c"
#include "unsource.h"
#undef  BASE_SHORT

#define BASE_UCHAR
#include "source.h"
#include "init_source.c"
#include "unsource.h"
#undef  BASE_UCHAR

#define BASE_CHAR
#include "source.h"
#include "init_source.c"
#include "unsource.h"
#undef  BASE_CHAR
