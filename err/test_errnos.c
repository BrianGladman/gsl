#include <config.h>

#include <stdio.h>
#include <stdlib.h>

#include <gsl_errno.h>
#include <gsl_test.h>

#define CHECK(x) errors[n].number = x ; errors[n].name = #x ; n++ ;
#define MAX_ERRS 32

int verbose = 0 ;

int
main (void)
{
  int i, j, n = 0 ;

  struct { 
    int number; 
    const char * name; 
  } errors[MAX_ERRS] ;

  CHECK(GSL_EDOM);
  CHECK(GSL_ERANGE);
  CHECK(GSL_EFAULT);
  CHECK(GSL_EINVAL);
  CHECK(GSL_EFAILED);
  CHECK(GSL_EFACTOR);
  CHECK(GSL_ESANITY);
  CHECK(GSL_ENOMEM);
  CHECK(GSL_EBADFUNC);
  CHECK(GSL_ERUNAWAY);
  CHECK(GSL_EMAXITER);
  CHECK(GSL_EZERODIV);
  CHECK(GSL_ETOL);
  CHECK(GSL_EUNDRFLW);
  CHECK(GSL_EOVRFLW);
  CHECK(GSL_ELOSS);

  for (i = 0 ; i < n ; i++) 
    {
      if (verbose) printf ("%s = %d\n", errors[i].name, errors[i].number) ;
    }

  for (i = 0; i < n; i++)
    {
      int status = 0;
      for (j = 0; j < n; j++)
	{
	  if (j != i)
	    status |= (errors[i].number == errors[j].number);
	}

      gsl_test (status, "%s is distinct from other error values",
		errors[i].name);
    }
  
  return gsl_test_summary ();
}

