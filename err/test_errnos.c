#include <config.h>

#include <stdio.h>
#include <stdlib.h>

#include <gsl_errno.h>
#include <gsl_test.h>

#define CHECK(x) {x,#x}
#define MAX_ERRS 16

int verbose = 0 ;

int
main (void)
{

  struct { 
    int number; 
    const char * name; 
  } errors[MAX_ERRS] = {
    CHECK(GSL_EDOM),
    CHECK(GSL_ERANGE),
    CHECK(GSL_EFAULT),
    CHECK(GSL_EINVAL),
    CHECK(GSL_EFAILED),
    CHECK(GSL_EFACTOR),
    CHECK(GSL_ESANITY),
    CHECK(GSL_ENOMEM),
    CHECK(GSL_EBADFUNC),
    CHECK(GSL_ERUNAWAY),
    CHECK(GSL_ETIMEOUT),
    CHECK(GSL_EZERODIV),
    CHECK(GSL_ETOL),
    {-1, "end"}
  } ;
  
  int i, j, n ;
  int status ;

  for (n = MAX_ERRS - 1 ; n > 0 ; n--) 
    {
      if (errors[n].number == -1) break ;
    }

  for (i = 0 ; i < n ; i++) 
    {
      if (verbose) printf ("%s = %d\n", errors[i].name, errors[i].number) ;
    }
  
  for (i = 0 ; i < n ; i++) 
    {
      status = (errors[i].number == 0)  ;
      gsl_test (status, "%s is non-zero", errors[i].name) ;

      for (j = 0 ; j < i ; j++) 
	{
	  status = (errors[i].number == errors[j].number) ;
	  gsl_test (status, "%s is distinct from %s", 
		    errors[j].name, errors[i].name) ;
	}
    }

  return gsl_test_summary ();

}

