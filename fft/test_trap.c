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
#include <gsl_fft_test.h>

#include <autotest.h>
#include <compare.h>

int verbose = 0;

unsigned int tests = 0;
unsigned int passed = 0;
unsigned int failed = 0;

void my_error_handler (const char *reason, const char *file, int line);

int
main (void)
{
  int status ;

  double real_x ;
  complex complex_x ;

  double * real_data = &real_x ;
  complex * complex_data = &complex_x  ; 

  gsl_fft_complex_wavetable complex_wavetable;
  gsl_fft_real_wavetable real_wavetable;
  gsl_fft_halfcomplex_wavetable halfcomplex_wavetable;

  gsl_set_error_handler (&my_error_handler);

  /* n = 0 in alloc */

  msg_checking ("trap for n = 0 in gsl_fft_complex_wavetable_alloc ");
  status = gsl_fft_complex_wavetable_alloc (0, &complex_wavetable);
  msg_result_status (!status);

  msg_checking ("trap for n = 0 in gsl_fft_real_wavetable_alloc ");
  status = gsl_fft_real_wavetable_alloc (0, &real_wavetable);
  msg_result_status (!status);

  msg_checking ("trap for n = 0 in gsl_fft_halfcomplex_wavetable_alloc ");
  status = gsl_fft_halfcomplex_wavetable_alloc (0, &halfcomplex_wavetable);
  msg_result_status (!status);


  /* n = 0 in wavetable_init */

  msg_checking ("trap for n = 0 in gsl_fft_complex_wavetable_init");
  status = gsl_fft_complex_init (0, &complex_wavetable);
  msg_result_status (!status);

  msg_checking ("trap for n = 0 in gsl_fft_real_wavetable_init");
  status = gsl_fft_real_init (0, &real_wavetable);
  msg_result_status (!status);

  msg_checking ("trap for n = 0 in gsl_fft_halfcomplex_wavetable_init");
  status = gsl_fft_halfcomplex_init (0, &halfcomplex_wavetable);
  msg_result_status (!status);


  /* n = 0 in generate_wavetable */

  msg_checking ("trap for n = 0 in gsl_fft_complex_generate_wavetable");
  status = gsl_fft_complex_generate_wavetable (0, &complex_wavetable);
  msg_result_status (!status);

  msg_checking ("trap for n = 0 in gsl_fft_real_generate_wavetable");
  status = gsl_fft_real_generate_wavetable (0, &real_wavetable);
  msg_result_status (!status);

  msg_checking ("trap for n = 0 in gsl_fft_halfcomplex_generate_wavetable");
  status = gsl_fft_halfcomplex_generate_wavetable (0, &halfcomplex_wavetable);
  msg_result_status (!status);


  /* n = 0 in fft forward */

  msg_checking ("trap for n = 0 in gsl_fft_complex_forward");
  status = gsl_fft_complex_forward (complex_data, 0, &complex_wavetable);
  msg_result_status (!status);

  msg_checking ("trap for n = 0 in gsl_fft_real");
  status = gsl_fft_real (real_data, 0, &real_wavetable);
  msg_result_status (!status);

  msg_checking ("trap for n = 0 in gsl_fft_halfcomplex");
  status = gsl_fft_halfcomplex (real_data, 0, &halfcomplex_wavetable);
  msg_result_status (!status);

  msg_checking ("trap for n = 0 in gsl_fft_complex_radix2_forward");
  status = gsl_fft_complex_radix2_forward (complex_data, 0);
  msg_result_status (!status);


  /* n = 0 in fft backward */

  msg_checking ("trap for n = 0 in gsl_fft_complex_backward");
  status = gsl_fft_complex_backward (complex_data, 0, &complex_wavetable);
  msg_result_status (!status);

  msg_checking ("trap for n = 0 in gsl_fft_complex_radix2_backward");
  status = gsl_fft_complex_radix2_backward (complex_data, 0);
  msg_result_status (!status);

  /* n = 0 in fft inverse */

  msg_checking ("trap for n = 0 in gsl_fft_complex_inverse");
  status = gsl_fft_complex_inverse (complex_data, 0, &complex_wavetable);
  msg_result_status (!status);

  msg_checking ("trap for n = 0 in gsl_fft_complex_radix2_inverse");
  status = gsl_fft_complex_radix2_inverse (complex_data, 0);
  msg_result_status (!status);


  /* n != 2^k in power of 2 routines */

  msg_checking ("trap for n != 2^k in gsl_fft_complex_radix2_forward");
  status = gsl_fft_complex_radix2_forward (complex_data, 17);
  msg_result_status (!status);

  msg_checking ("trap for n != 2^k in gsl_fft_complex_radix2_backward");
  status = gsl_fft_complex_radix2_backward (complex_data, 17);
  msg_result_status (!status);

  msg_checking ("trap for n != 2^k in gsl_fft_complex_radix2_inverse");
  status = gsl_fft_complex_radix2_inverse (complex_data, 17);
  msg_result_status (!status);

  /* n != wavetable.n in mixed radix routines */

  msg_checking ("trap for n != nw in gsl_fft_complex_forward");
  complex_wavetable.n = 3;
  status = gsl_fft_complex_forward (complex_data, 4, &complex_wavetable);
  msg_result_status (!status);

  msg_checking ("trap for n != nw in gsl_fft_complex_backward");
  complex_wavetable.n = 3;
  status = gsl_fft_complex_backward (complex_data, 4, &complex_wavetable);
  msg_result_status (!status);

  msg_checking ("trap for n != nw in gsl_fft_complex_inverse");
  complex_wavetable.n = 3;
  status = gsl_fft_complex_inverse (complex_data, 4, &complex_wavetable);
  msg_result_status (!status);

  msg_checking ("trap for n != nw in gsl_fft_real");
  real_wavetable.n = 3;
  status = gsl_fft_real (real_data, 4, &real_wavetable);
  msg_result_status (!status);

  msg_checking ("trap for n != nw in gsl_fft_halfcomplex");
  halfcomplex_wavetable.n = 3;
  status = gsl_fft_halfcomplex (real_data, 4, &halfcomplex_wavetable);
  msg_result_status (!status);

  status = msg_summary (tests, passed, failed);
  
  return status;
}


void
my_error_handler (const char *reason, const char *file, int line)
{
  printf ("(caught)\n[%s:%d: %s]\n", file, line, reason);
}
