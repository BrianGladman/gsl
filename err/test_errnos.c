#include <config.h>

#include <stdio.h>
#include <stdlib.h>

#include <gsl_errno.h>

#define CHECK(x) {x,#x}
#define MAX_ERRS 16

main () {

  struct { 
    int number; 
    const char * name; 
  } errors[MAX_ERRS] = {  CHECK(GSL_EDOM),
			  CHECK(GSL_ERANGE),
			  CHECK(GSL_EFAULT),
			  CHECK(GSL_EINVAL),
			  CHECK(GSL_EFAILED),
			  CHECK(GSL_EFACTOR),
			  CHECK(GSL_ESANITY),
			  CHECK(GSL_ENOMEM),
			  -1 } ;
  
  int errs, i, j ;
  int failed = 0 ;

  for (errs = MAX_ERRS - 1 ; errs > 0 ; errs--) 
    {
      if (errors[errs].number == -1) break ;
    }

  for (i = 0 ; i < errs ; i++) 
    {
      printf ("%s = %d\n", errors[i].name, errors[i].number) ;
    }
  
  for (i = 0 ; i < errs ; i++) 
    {
      if (errors[i].number == 0) 
	{
	  printf ("Error: %s is 0\n", errors[i].name) ;
	  failed++ ;
	}
      for (j = 0 ; j < i ; j++) 
	{
	  if (errors[i].number == errors[j].number) {
	    printf ("Error: %s is ambiguous with %s (both = %d)\n", 
		    errors[j].name, errors[i].name, errors[i].number) ;
	    failed++ ;
	  }
	}
    }

  if (failed) 
    {
      exit(EXIT_FAILURE) ;
    } 
  else 
    {
      printf("Error codes are unique and non-zero\n"); 
      exit (EXIT_SUCCESS) ;
    }

}

