#include <config.h>

#include <stdio.h>
#include <stdlib.h>

#if HAVE_VPRINTF
#if __STDC__
#include <stdarg.h>
#else
#include <varargs.h>
#endif
#endif

#include <gsl_test.h>

static unsigned int tests = 0 ;
static unsigned int passed = 0 ;
static unsigned int failed = 0 ; 

static unsigned int verbose = 1 ;

void
gsl_test (int status, const char *test_description, ...)
{
  
  tests++;
  
  if (status == 0)
    {
      passed++;
      if (verbose) printf ("PASS: ");
    }
  else
    {
      failed++;
      if (verbose) printf ("FAIL: ");
    }

  if (verbose) 
    {

#ifdef HAVE_VPRINTF
      va_list ap ;

#if __STDC__
      va_start (ap, test_description) ;
#else
      va_start (ap) ;
#endif
      vprintf (test_description, ap);
      va_end (ap) ;
#endif

      putchar ('\n') ;
      fflush (stdout);
    }
}

void 
gsl_test_verbose (int v)
{
  verbose = v ;
}

int
gsl_test_summary ()
{

  if (verbose && 0) /* FIXME: turned it off, this annoys me */
    printf ("%d tests, passed %d, failed %d.\n", tests, passed, failed);

  if (failed != 0)
    {

      if (verbose && 0)  /* FIXME: turned it off, this annoys me */
	{
	  printf ("%d TEST%s FAILED.\n", failed, failed == 1 ? "" : "S" );
	}
      return EXIT_FAILURE ;
    }

  if (tests != passed + failed)
    {
      if (verbose) printf ("TEST RESULTS DO NOT ADD UP %d != %d + %d\n",
			   tests, passed, failed);
      return EXIT_FAILURE ;
    }

  if (passed == tests)  
    {
      if (verbose && 0) /* FIXME: turned it off, this annoys me */
	printf ("All tests passed successfully\n");
      return EXIT_SUCCESS ;
    }
  
  return EXIT_FAILURE;
}
