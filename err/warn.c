#include <config.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>

#include <gsl_errno.h>
#include <gsl_message.h>

#ifdef GSL_THREAD_SAFE
void
gsl_warning (const char * reason, const char * file, int line, int gsl_errno)
{
  gsl_errno = 0;		/* stop complaints about unused variables */
  gsl_stream_printf ("WARNING", file, line, reason);
}
#else /* GSL_THREAD_SAFE */
int gsl_warnings_off = 0;

void
gsl_warning (const char * reason, const char * file, int line, int gsl_errno)
{
  if (!gsl_warnings_off)
    {
      const char * error_string = gsl_strerror(gsl_errno) ;
      gsl_errno = 0;		/* stop complaints about unused variables */
      gsl_stream_printf ("WARNING", file, line, reason);
      gsl_stream_printf ("(ERRNO)", file, line, error_string);
    }
}
#endif /* GSL_THREAD_SAFE */


