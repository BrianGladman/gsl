#ifndef GSL_TEST_H
#define GSL_TEST_H

void
  gsl_test (int status, const char *test_description, ...);

void
gsl_test_rel (double result, double expected, double relative_error,
	      const char *test_description, ...) ;

void
gsl_test_abs (double result, double expected, double relative_error,
	      const char *test_description, ...) ;

void
gsl_test_int (int result, int expected, const char *test_description, ...) ;

void
gsl_test_str (const char * result, const char * expected, 
	      const char *test_description, ...) ;

void
  gsl_test_verbose (int verbose) ;

int
  gsl_test_summary (void) ;


#endif /* GSL_TEST_H */
