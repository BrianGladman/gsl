#ifndef __GSL_TEST_H__
#define __GSL_TEST_H__

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

void
  gsl_test (int status, const char *test_description, ...);

void
gsl_test_rel (double result, double expected, double relative_error,
	      const char *test_description, ...) ;

void
gsl_test_abs (double result, double expected, double absolute_error,
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


__END_DECLS

#endif /* __GSL_TEST_H__ */
