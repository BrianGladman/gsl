#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_heapsort.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>

double urand (void);

#include "test_heapsort.c"

#define BASE_LONG_DOUBLE
#include "templates_on.h"
#include "test_source.c"
#include "templates_off.h"
#undef  BASE_LONG_DOUBLE

#define BASE_DOUBLE
#include "templates_on.h"
#include "test_source.c"
#include "templates_off.h"
#undef  BASE_DOUBLE

#define BASE_FLOAT
#include "templates_on.h"
#include "test_source.c"
#include "templates_off.h"
#undef  BASE_FLOAT

#define BASE_ULONG
#include "templates_on.h"
#include "test_source.c"
#include "templates_off.h"
#undef  BASE_ULONG

#define BASE_LONG
#include "templates_on.h"
#include "test_source.c"
#include "templates_off.h"
#undef  BASE_LONG

#define BASE_UINT
#include "templates_on.h"
#include "test_source.c"
#include "templates_off.h"
#undef  BASE_UINT

#define BASE_INT
#include "templates_on.h"
#include "test_source.c"
#include "templates_off.h"
#undef  BASE_INT

#define BASE_USHORT
#include "templates_on.h"
#include "test_source.c"
#include "templates_off.h"
#undef  BASE_USHORT

#define BASE_SHORT
#include "templates_on.h"
#include "test_source.c"
#include "templates_off.h"
#undef  BASE_SHORT

#define BASE_UCHAR
#include "templates_on.h"
#include "test_source.c"
#include "templates_off.h"
#undef  BASE_UCHAR

#define BASE_CHAR
#include "templates_on.h"
#include "test_source.c"
#include "templates_off.h"
#undef  BASE_CHAR

int
main (void)
{
  size_t i, s;

  /* Test for lengths of 1 ... 31, then 32, 64, 128, 256, ... */

  for (i = 1; i < 1024; i = (i < 32) ? i + 1 : 2 * i)
    test_heapsort (i);

  for (i = 1; i < 1024; i = (i < 32) ? i + 1 : 2 * i)
    {
      for (s = 1; s < 4; s++)
	{
	  test_sort_vector (i, s);
	  test_sort_vector_float (i, s);
	  test_sort_vector_long_double (i, s);
	  test_sort_vector_ulong (i, s);
	  test_sort_vector_long (i, s);
	  test_sort_vector_uint (i, s);
	  test_sort_vector_int (i, s);
	  test_sort_vector_ushort (i, s);
	  test_sort_vector_short (i, s);
	  test_sort_vector_uchar (i, s);
	  test_sort_vector_char (i, s);
	}
    }

  return gsl_test_summary ();
}

double 
urand (void)
{
  static unsigned long int x = 1;
  x = (1103515245 * x + 12345) & 0x7fffffffUL;
  return x / 2147483648.0;
}
