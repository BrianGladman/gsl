#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>

#include <gsl_errno.h>

int gsl_errno;

typedef void handler (const char *reason, const char *file, int line);

handler *gsl_error_handler = NULL;

handler *
gsl_error_set_handler (handler * new_handler)
{
  handler *previous_handler = gsl_error_handler;
  gsl_error_handler = new_handler;
  return previous_handler;
}

void
gsl_error (const char *reason, const char *file, int line)
{
  if (gsl_error_handler)
    {
      (*gsl_error_handler) (reason, file, line);
      return;
    }

  fprintf (stderr, "%s:%d: error: %s\n", file, line, reason);
  abort ();
}
