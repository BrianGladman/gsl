#include <config.h>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include <gsl_complex.h>
#include <gsl_errno.h>
#include <gsl_dft_complex.h>
#include <gsl_fft_complex.h>
#include <gsl_fft_real.h>
#include <gsl_fft_halfcomplex.h>
#include <gsl_test.h>

#include "fft_signals.h"

#include "compare.h"

int verbose = 0;

void usage (void);
void my_error_handler (const char *reason, const char *file,
		       int line, int err);

/* Usage: test_trap
   Exercise the error handling in the fft routines.
   The exit status indicates success or failure. */

int
main (void)
{
  int status ;

  double real_x ;
  gsl_complex complex_x ;

  double * real_data = &real_x ;
  gsl_complex * complex_data = &complex_x  ; 

  gsl_fft_complex_wavetable * cw;
  gsl_fft_real_wavetable * rw;
  gsl_fft_halfcomplex_wavetable * hcw;

  gsl_set_error_handler (&my_error_handler);

  /* n = 0 in alloc */

  cw = gsl_fft_complex_wavetable_alloc (0);
  gsl_test (cw != 0, "trap for n = 0 in gsl_fft_complex_wavetable_alloc");

  rw = gsl_fft_real_wavetable_alloc (0);
  gsl_test (rw != 0, "trap for n = 0 in gsl_fft_real_wavetable_alloc" );

  hcw = gsl_fft_halfcomplex_wavetable_alloc (0);
  gsl_test (hcw != 0, "trap for n = 0 in gsl_fft_halfcomplex_wavetable_alloc");

  /* n = 0 in wavetable_init */

  status = gsl_fft_complex_init (0, cw);
  gsl_test (!status, "trap for n = 0 in gsl_fft_complex_wavetable_init");

  status = gsl_fft_real_init (0, rw);
  gsl_test (!status, "trap for n = 0 in gsl_fft_real_wavetable_init");

  status = gsl_fft_halfcomplex_init (0, hcw);
  gsl_test (!status, "trap for n = 0 in gsl_fft_halfcomplex_wavetable_init");


  /* n = 0 in generate_wavetable */

  status = gsl_fft_complex_generate_wavetable (0, cw);
  gsl_test (!status, "trap for n = 0 in gsl_fft_complex_generate_wavetable");

  status = gsl_fft_real_generate_wavetable (0, rw);
  gsl_test (!status, "trap for n = 0 in gsl_fft_real_generate_wavetable");

  status = gsl_fft_halfcomplex_generate_wavetable (0, hcw);
  gsl_test (!status, "trap for n = 0 in gsl_fft_halfcomplex_generate_wavetable");

  cw = gsl_fft_complex_wavetable_alloc (10);
  hcw = gsl_fft_halfcomplex_wavetable_alloc (10);
  rw = gsl_fft_real_wavetable_alloc (10);

  /* n = 0 in fft forward */

  status = gsl_fft_complex_forward (complex_data, 0, cw);
  gsl_test (!status, "trap for n = 0 in gsl_fft_complex_forward");

  status = gsl_fft_real (real_data, 0, rw);
  gsl_test (!status, "trap for n = 0 in gsl_fft_real");

  status = gsl_fft_halfcomplex (real_data, 0, hcw);
  gsl_test (!status, "trap for n = 0 in gsl_fft_halfcomplex");

  status = gsl_fft_complex_radix2_forward (complex_data, 0);
  gsl_test (!status, "trap for n = 0 in gsl_fft_complex_radix2_forward");

  /* n = 0 in fft backward */

  status = gsl_fft_complex_backward (complex_data, 0, cw);
  gsl_test (!status, "trap for n = 0 in gsl_fft_complex_backward");

  status = gsl_fft_complex_radix2_backward (complex_data, 0);
  gsl_test (!status, "trap for n = 0 in gsl_fft_complex_radix2_backward");

  /* n = 0 in fft inverse */

  status = gsl_fft_complex_inverse (complex_data, 0, cw);
  gsl_test (!status, "trap for n = 0 in gsl_fft_complex_inverse");

  status = gsl_fft_complex_radix2_inverse (complex_data, 0);
  gsl_test (!status, "trap for n = 0 in gsl_fft_complex_radix2_inverse");

  /* n != 2^k in power of 2 routines */

  status = gsl_fft_complex_radix2_forward (complex_data, 17);
  gsl_test (!status, "trap for n != 2^k in gsl_fft_complex_radix2_forward");

  status = gsl_fft_complex_radix2_backward (complex_data, 17);
  gsl_test (!status, "trap for n != 2^k in gsl_fft_complex_radix2_backward");

  status = gsl_fft_complex_radix2_inverse (complex_data, 17);
  gsl_test (!status, "trap for n != 2^k in gsl_fft_complex_radix2_inverse");

  /* n != wavetable.n in mixed radix routines */

  cw->n = 3;
  status = gsl_fft_complex_forward (complex_data, 4, cw);
  gsl_test (!status, "trap for n != nw in gsl_fft_complex_forward");

  cw->n = 3;
  status = gsl_fft_complex_backward (complex_data, 4, cw);
  gsl_test (!status, "trap for n != nw in gsl_fft_complex_backward");

  cw->n = 3;
  status = gsl_fft_complex_inverse (complex_data, 4, cw);
  gsl_test (!status, "trap for n != nw in gsl_fft_complex_inverse");

  rw->n = 3;
  status = gsl_fft_real (real_data, 4, rw);
  gsl_test (!status, "trap for n != nw in gsl_fft_real");

  hcw->n = 3;
  status = gsl_fft_halfcomplex (real_data, 4, hcw);
  gsl_test (!status, "trap for n != nw in gsl_fft_halfcomplex");

  return gsl_test_summary ();

}


void
my_error_handler (const char *reason, const char *file, int line, int err)
{
  if (0) printf ("(caught [%s:%d: %s (%d)])\n", file, line, reason, err) ;
}
