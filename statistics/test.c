#include <config.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_test.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_statistics.h>

/* Test program for mean.c.  JimDavies 7.96 */

#define BASE_LONG_DOUBLE
#include "templates_on.h"
#include "test_float_source.c"
#include "templates_off.h"
#undef  BASE_LONG_DOUBLE

#define BASE_DOUBLE
#include "templates_on.h"
#include "test_float_source.c"
#include "templates_off.h"
#undef  BASE_DOUBLE

#define BASE_FLOAT
#include "templates_on.h"
#include "test_float_source.c"
#include "templates_off.h"
#undef  BASE_FLOAT

#define BASE_ULONG
#include "templates_on.h"
#include "test_int_source.c"
#include "templates_off.h"
#undef  BASE_ULONG

#define BASE_LONG
#include "templates_on.h"
#include "test_int_source.c"
#include "templates_off.h"
#undef  BASE_LONG

#define BASE_UINT
#include "templates_on.h"
#include "test_int_source.c"
#include "templates_off.h"
#undef  BASE_UINT

#define BASE_INT
#include "templates_on.h"
#include "test_int_source.c"
#include "templates_off.h"
#undef  BASE_INT

#define BASE_USHORT
#include "templates_on.h"
#include "test_int_source.c"
#include "templates_off.h"
#undef  BASE_USHORT

#define BASE_SHORT
#include "templates_on.h"
#include "test_int_source.c"
#include "templates_off.h"
#undef  BASE_SHORT

#define BASE_UCHAR
#include "templates_on.h"
#include "test_int_source.c"
#include "templates_off.h"
#undef  BASE_UCHAR

#define BASE_CHAR
#include "templates_on.h"
#include "test_int_source.c"
#include "templates_off.h"
#undef  BASE_CHAR


int
main (void)
{
  size_t s1, s2;

  for (s1 = 1; s1 < 4 ; s1++) 
    {
      for (s2 = 1; s2 < 4 ; s2++) 
        {
          test_func (s1,s2);
          test_float_func (s1,s2);
          test_long_double_func (s1,s2);
          
          test_ulong_func (s1,s2);
          test_long_func (s1,s2);
          test_uint_func (s1,s2);
          test_int_func (s1,s2);
          test_ushort_func (s1,s2);
          test_short_func (s1,s2);
          test_uchar_func (s1,s2);
          test_char_func (s1,s2);
        }
    }

  return gsl_test_summary ();
}

