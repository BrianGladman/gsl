#include <config.h>

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include <gsl_complex.h>
#include <gsl_errno.h>
#include <gsl_dft_complex.h>
#include <gsl_fft_complex.h>
#include <gsl_fft_complex_float.h>
#include <gsl_fft_real.h>
#include <gsl_fft_real_float.h>
#include <gsl_fft_halfcomplex.h>
#include <gsl_fft_halfcomplex_float.h>
#include <gsl_test.h>

void my_error_handler (const char *reason, const char *file,
		       int line, int err);

#include "complex_internal.h"

/* Usage: test [n]
   Exercise the fft routines for length n. By default n runs from 1 to 100.
   The exit status indicates success or failure. */

#define BASE_DOUBLE
#include "templates_on.h"
#include "test_complex.c"
#include "test_real.c"
#include "test_trap.c"
#include "templates_off.h"
#undef  BASE_DOUBLE

#define BASE_FLOAT
#include "templates_on.h"
#include "test_complex.c"
#include "test_real.c"
#include "test_trap.c"
#include "templates_off.h"
#undef  BASE_FLOAT

int
main (int argc, char *argv[])
{
  size_t i;
  size_t start = 1, end = 99;
  size_t stride ;
  size_t n = 0;

  if (argc == 2) 
    n = strtol (argv[1], NULL, 0);
 
  if (n)
    {
      start = n ;
      end = n ;
    }


  for (i = 1 ; i <= end ; i *= 2) 
    {
      if (i >= start) 
	{
	  for (stride = 1 ; stride < 4 ; stride++)
	    {
	      test_complex_bitreverse_order (stride, i) ;
	      test_complex_radix2 (stride, i) ;
	      test_real_bitreverse_order (stride, i) ;
	      test_real_radix2 (stride, i) ;
	    }
	}
    }

  for (i = start ; i <= end ; i++) 
    {
      for (stride = 1 ; stride < 4 ; stride++)
	{
	  test_complex_func (stride, i) ;
	  test_complex_float_func (stride, i) ;
	  test_real_func (stride, i) ;
	  test_real_float_func (stride, i) ;
	}
    }

  gsl_set_error_handler (&my_error_handler);
  test_trap () ;
  test_float_trap () ;

  return gsl_test_summary ();
}


void
my_error_handler (const char *reason, const char *file, int line, int err)
{
  if (0) printf ("(caught [%s:%d: %s (%d)])\n", file, line, reason, err) ;
}
