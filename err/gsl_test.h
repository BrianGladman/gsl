#ifndef GSL_TEST_H
#define GSL_TEST_H

#include <config.h>

#ifdef HAVE_VPRINTF
void
  gsl_test (int status, const char *test_description, ...);
#endif 

void
  gsl_test_verbose (int verbose) ;

int
  gsl_test_summary (void) ;


#endif /* GSL_TEST_H */
