#include <config.h>

#include <stdio.h>
#include <stdlib.h>

#include <autotest.h>

extern unsigned int tests;
extern unsigned int passed;
extern unsigned int failed;

void
msg_checking (const char *test_description)
{
  printf ("checking %s... ", test_description);
  fflush (stdout);
}

void
msg_checking_params (const char *params_description,
		     const char *test_description)
{
  printf ("%s checking %s... ", params_description, test_description);
  fflush (stdout);
}


void
msg_result_status (int status)
{
  tests++;

  if (status == 0)
    {
      passed++;
      msg_result ("ok");
    }
  else
    {
      failed++;
      msg_error ("failed");
    }
}

void
msg_result (const char *result_description)
{
  printf ("%s\n", result_description);
}

void
msg_error (const char *error_description)
{
  printf ("%s\n", error_description);
}

int
msg_summary (unsigned int total_tests,
	     unsigned int total_passed,
	     unsigned int total_failed)
{

  printf ("%d tests. Passed %d. Failed %d.\n",
	  total_tests, total_passed, total_failed);

  if (total_failed != 0)
    {
      printf ("%d TEST(S) FAILED.\n", total_failed);
      return EXIT_FAILURE ;
    }

  if (total_tests != total_passed + total_failed)
    {
      printf ("Test results do not add up %d != %d + %d\n",
	      total_tests, total_passed, total_failed);
      return EXIT_FAILURE ;
    }

  if (total_passed == total_tests)
    {
      printf ("All tests passed successfully\n");
      return EXIT_SUCCESS ;
    }

  return EXIT_FAILURE;
}
