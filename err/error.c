#include <config.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>

#include <gsl_errno.h>
#include <gsl_message.h>

#ifdef GSL_THREAD_SAFE
void
gsl_error (const char * reason, const char * file, int line, int gsl_errno)
  gsl_stream_printf ("ERROR", file, line, reason);
  abort ();
}
#else /* GSL_THREAD_SAFE */
gsl_error_handler_t * gsl_error_handler = NULL;

void
gsl_error (const char * reason, const char * file, int line, int gsl_errno)
{
  if (gsl_error_handler) 
    {
      (*gsl_error_handler) (reason, file, line, gsl_errno);
      return ;
    }

  gsl_stream_printf ("ERROR", file, line, reason);
  abort ();
}

gsl_error_handler_t *
gsl_set_error_handler (gsl_error_handler_t * new_handler)
{
  gsl_error_handler_t * previous_handler = gsl_error_handler;
  gsl_error_handler = new_handler;
  return previous_handler;
}
#endif /* GSL_THREAD_SAFE */

