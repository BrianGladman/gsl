#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>

#include <gsl_errno.h>

int gsl_errno;

gsl_errhandler_t *gsl_error_handler = NULL;

gsl_errhandler_t *
gsl_set_error_handler (gsl_errhandler_t * new_handler)
{
  gsl_errhandler_t *previous_handler = gsl_error_handler;
  gsl_error_handler = new_handler;
  return previous_handler;
}

void 
gsl_no_error_handler (const char * reason, const char * file, int line) {
  reason = 0 ; /* stop complaints about unused parameters */
  file = 0 ;
  line = 0 ;
}

void
gsl_error (const char *reason, const char *file, int line)
{
  if (gsl_error_handler)
    {
      (*gsl_error_handler) (reason, file, line);
      return;
    }

  fprintf (stderr, "\n%s:%d: error: %s\n", file, line, reason);
  abort ();
}
