#include <config.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>

#include <gsl_errno.h>
#include <gsl_message.h>

static void 
gsl_error_stream_printf(const char * label, const char * file, 
			int line, const char * reason) ;

int gsl_warnings_off = 0 ;
unsigned int gsl_message_mask = GSL_MESSAGE_MASK;


gsl_errhandler_t *gsl_error_handler = NULL;

gsl_errhandler_t *
gsl_set_error_handler (gsl_errhandler_t * new_handler)
{
  gsl_errhandler_t *previous_handler = gsl_error_handler;
  gsl_error_handler = new_handler;
  return previous_handler;
}

void 
gsl_empty_error_handler (const char * reason, const char * file, int line,
			 int gsl_errno) 
{
  reason = 0 ; /* stop complaints about unused parameters */
  file = 0 ;
  line = 0 ;
  gsl_errno = 0 ;
}

static void
gsl_error_stream_printf(const char * label, const char * file, 
			int line, const char * reason)
{
  fprintf(stderr, "\ngsl: %s:%d: %s: %s\n", file, line, label, reason);
}

void
gsl_warning(const char * reason, const char * file, int line)
{
  if(!gsl_warnings_off) 
    gsl_error_stream_printf("WARNING", file, line, reason);
}

void
gsl_error (const char *reason, const char *file, int line, int gsl_errno)
{
  if (gsl_error_handler)
    {
      (*gsl_error_handler) (reason, file, line, gsl_errno);
      return;
    }

  gsl_error_stream_printf("ERROR", file, line, reason);
  abort ();
}

void
gsl_message(const char * reason, const char * file,
	    int line, unsigned int mask)
{
  if(mask & gsl_message_mask) {
    gsl_error_stream_printf("MESSAGE", file, line, reason);
  }
}
