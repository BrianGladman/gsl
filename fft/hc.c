#include <config.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_complex.h>

#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_fft_halfcomplex_float.h>

#define BASE_DOUBLE
#include "templates_on.h"
#include "hc_pass.h"
#include "hc_init.c"
#include "hc_main.c"
#include "hc_pass_2.c"
#include "hc_pass_3.c"
#include "hc_pass_4.c"
#include "hc_pass_5.c"
#include "hc_pass_n.c"
#include "hc_radix2.c"
#include "hc_unpack.c"
#include "templates_off.h"
#undef  BASE_DOUBLE

#define BASE_FLOAT
#include "templates_on.h"
#include "hc_pass.h"
#include "hc_init.c"
#include "hc_main.c"
#include "hc_pass_2.c"
#include "hc_pass_3.c"
#include "hc_pass_4.c"
#include "hc_pass_5.c"
#include "hc_pass_n.c"
#include "hc_radix2.c"
#include "hc_unpack.c"
#include "templates_off.h"
#undef  BASE_FLOAT
