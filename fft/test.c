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
#include <gsl_fft_halfcomplex.h>
#include <gsl_test.h>

#include "complex_internal.h"

/* Usage: test [n]
   Exercise the fft routines for length n. By default n runs from 1 to 100.
   The exit status indicates success or failure. */

#define BASE_DOUBLE
#include "templates_on.h"
#include "test_complex.c"
#include "templates_off.h"
#undef  BASE_DOUBLE

#define BASE_FLOAT
#include "templates_on.h"
#include "test_complex.c"
#include "templates_off.h"
#undef  BASE_FLOAT

int
main (int argc, char *argv[])
{
  size_t i, start = 1, end = 99;
  size_t n = 0;

  if (argc == 2) 
    n = strtol (argv[1], NULL, 0);
 
  if (n)
    {
      start = n ;
      end = n ;
    }

  for (i = start ; i < end + 1; i++) 
    {
      test_complex_func (i, 1) ;
      test_complex_float_func (i, 1) ;
      test_complex_func (i, 3) ;
      test_complex_float_func (i, 3) ;
    }


  return gsl_test_summary ();
}

