#include <config.h>
#include <stdio.h>

#undef GSL_RANGE_CHECK_OFF	/* we want range checking for this test */
#define GSL_WARNINGS_OFF

#include <gsl_vector.h>
#include <gsl_test.h>

#define N 10000

void my_error_handler (const char *reason, const char *file,
		       int line, int err);

int status = 0;

int 
main (void)
{
  gsl_vector *v = gsl_vector_alloc (N);
  gsl_vector_int *vi = gsl_vector_int_alloc (N);
  gsl_vector_float *vf = gsl_vector_float_alloc (N);
  gsl_vector_complex *vc = gsl_vector_complex_alloc (N);
  size_t j = 0;
  double x;
  gsl_complex z1;
  gsl_complex z =
  {
    {1.2, 3.4}};

  gsl_warnings_off = 1;

  gsl_set_error_handler (&my_error_handler);

  status = 0;
  gsl_vector_set (v, j - 1, 1.2);
  gsl_test (!status, "gsl_vector_set traps index below lower array bound");

  status = 0;
  gsl_vector_set (v, N + 1, 1.2);
  gsl_test (!status, "gsl_vector_set traps index above upper array bound");

  status = 0;
  gsl_vector_set (v, N, 1.2);
  gsl_test (!status, "gsl_vector_set traps index at upper array bound");

  status = 0;
  x = gsl_vector_get (v, j - 1);
  gsl_test (!status, "gsl_vector_get traps index below lower array bound");
  gsl_test(x != 0,
	   "gsl_vector_get returns zero for index below lower array bound") ;

  status = 0;
  x = gsl_vector_get (v, N + 1);
  gsl_test (!status, "gsl_vector_get traps index above upper array bound");
  gsl_test(x != 0,
	   "gsl_vector_get returns zero for index above upper array bound") ;

  status = 0;
  x = gsl_vector_get (v, N);
  gsl_test (!status, "gsl_vector_get traps index at upper array bound");
  gsl_test(x != 0,
	   "gsl_vector_get returns zero for index at upper array bound") ;

  /* integer */

  status = 0;
  gsl_vector_int_set (vi, j - 1, 1);
  gsl_test (!status, "gsl_vector_int_set traps index below lower array bound");

  status = 0;
  gsl_vector_int_set (vi, N + 1, 1);
  gsl_test (!status, "gsl_vector_int_set traps index above upper array bound");

  status = 0;
  gsl_vector_int_set (vi, N, 1);
  gsl_test (!status, "gsl_vector_int_set traps index at upper array bound");

  status = 0;
  x = gsl_vector_int_get (vi, j - 1);
  gsl_test (!status, "gsl_vector_int_get traps index below lower array bound");
  gsl_test(x != 0,
	   "gsl_vector_int_get returns zero below lower array bound") ;

  status = 0;
  x = gsl_vector_int_get (vi, N + 1);
  gsl_test (!status, "gsl_vector_int_get traps index above upper array bound");
  gsl_test(x != 0,
	   "gsl_vector_int_get returns zero above upper array bound") ;

  status = 0;
  x = gsl_vector_int_get (vi, N);
  gsl_test (!status, "gsl_vector_int_get traps index at upper array bound");
  gsl_test(x != 0,
	   "gsl_vector_int_get returns zero for index at upper array bound") ;

  /* float */

  status = 0;
  gsl_vector_float_set (vf, j - 1, 1.2);
  gsl_test (!status,
	    "gsl_vector_float_set traps index below lower array bound");

  status = 0;
  gsl_vector_float_set (vf, N + 1, 1.2);
  gsl_test (!status,
	    "gsl_vector_float_set traps index above upper array bound");

  status = 0;
  gsl_vector_float_set (vf, N, 1.2);
  gsl_test (!status, "gsl_vector_float_set traps index at upper array bound");

  status = 0;
  x = gsl_vector_float_get (vf, j - 1);
  gsl_test (!status,
	    "gsl_vector_float_get traps index below lower array bound");
  gsl_test(x != 0,
	   "gsl_vector_float_get returns zero below lower array bound") ;

  status = 0;
  x = gsl_vector_float_get (vf, N + 1);
  gsl_test (!status,
	    "gsl_vector_float_get traps index above upper array bound");
  gsl_test(x != 0,
	   "gsl_vector_float_get returns zero above upper array bound") ;

  status = 0;
  x = gsl_vector_float_get (vf, N);
  gsl_test (!status, "gsl_vector_float_get traps index at upper array bound");
  gsl_test(x != 0,
	   "gsl_vector_float_get returns zero at upper array bound") ;

  /* complex */

  status = 0;
  gsl_vector_complex_set (vc, j - 1, z);
  gsl_test (!status,
	    "gsl_vector_complex_set traps index below lower array bound");

  status = 0;
  gsl_vector_complex_set (vc, N + 1, z);
  gsl_test (!status,
	    "gsl_vector_complex_set traps index above upper array bound");

  status = 0;
  gsl_vector_complex_set (vc, N, z);
  gsl_test (!status, "gsl_vector_complex_set traps index at upper array bound");

  status = 0;
  z1 = gsl_vector_complex_get (vc, j - 1);
  gsl_test (!status,
	    "gsl_vector_complex_get traps index below lower array bound");

  gsl_test(GSL_COMPLEX_REAL(z1) != 0,
	   "gsl_vector_complex_get returns zero real below lower array bound") ;
  gsl_test(GSL_COMPLEX_IMAG(z1) != 0,
	   "gsl_vector_complex_get returns zero imag below lower array bound") ;

  status = 0;
  z1 = gsl_vector_complex_get (vc, N + 1);
  gsl_test (!status,
	    "gsl_vector_complex_get traps index above upper array bound");
  gsl_test(GSL_COMPLEX_REAL(z1) != 0,
	   "gsl_vector_complex_get returns zero real above upper array bound") ;
  gsl_test(GSL_COMPLEX_IMAG(z1) != 0,
	   "gsl_vector_complex_get returns zero imag above upper array bound") ;

  status = 0;
  z1 = gsl_vector_complex_get (vc, N);
  gsl_test (!status, "gsl_vector_complex_get traps index at upper array bound");
  gsl_test(GSL_COMPLEX_REAL(z1) != 0,
	   "gsl_vector_complex_get returns zero real at upper array bound") ;
  gsl_test(GSL_COMPLEX_IMAG(z1) != 0,
	   "gsl_vector_complex_get returns zero imag at upper array bound") ;

  return gsl_test_summary ();
}

void
my_error_handler (const char *reason, const char *file, int line, int err)
{
  if (0)
    printf ("(caught [%s:%d: %s (%d)])\n", file, line, reason, err);
  status = 1;
}
