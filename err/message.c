#include <config.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_message.h>

unsigned int gsl_message_mask = GSL_MESSAGE_MASK;

void
gsl_message (const char * reason, const char * file, int line, 
	     unsigned int mask)
{
  if (mask & gsl_message_mask)
    {
      gsl_stream_printf ("MESSAGE", file, line, reason);
    }
}
