#include <config.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_message.h>

FILE * gsl_stream = NULL ;
gsl_stream_handler_t * gsl_stream_handler = NULL;

void
gsl_stream_printf (const char *label, const char *file, int line, 
		   const char *reason)
{
  if (gsl_stream == NULL)
    {
      gsl_stream = stderr;
    }
  if (gsl_stream_handler)
    {
      (*gsl_stream_handler) (label, file, line, reason);
      return;
    }
  fprintf (gsl_stream, "gsl: %s:%d: %s: %s\n", file, line, label, reason);

}

gsl_stream_handler_t *
gsl_set_stream_handler (gsl_stream_handler_t * new_handler)
{
  gsl_stream_handler_t * previous_handler = gsl_stream_handler;
  gsl_stream_handler = new_handler;
  return previous_handler;
}

FILE *
gsl_set_stream (FILE * new_stream)
{
  FILE * previous_stream;
  if (gsl_stream == NULL) {
    gsl_stream = stderr;
  }
  previous_stream = gsl_stream;
  gsl_stream = new_stream;
  return previous_stream;
}
