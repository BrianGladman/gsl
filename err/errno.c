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
gsl_empty_error_handler (const char * reason, const char * file, int line) {
  reason = 0 ; /* stop complaints about unused parameters */
  file = 0 ;
  line = 0 ;
}

void
gsl_message(const char * reason, const char * file, int line, unsigned int mask)
{
  if(mask & GSL_MESSAGE_MASK) {
    fprintf(stderr, "\ngsl: %s:%d: MESSAGE: %s\n", file, line, reason);
  }
}

void
gsl_warning(const char * reason, const char * file, int line)
{
  fprintf(stderr, "\ngsl: %s:%d: WARNING: %s\n", file, line, reason);
}

void
gsl_error (const char *reason, const char *file, int line)
{
  if (gsl_error_handler)
    {
      (*gsl_error_handler) (reason, file, line);
      return;
    }

  gsl_message(reason, file, line);
  abort ();
}
