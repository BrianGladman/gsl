#include <config.h>
#include <stdlib.h>
#include <gsl_ieee_utils.h>

static int little_endian_p (void) ;

static int 
little_endian_p (void) {
  /* Are we little or big endian?  From Harbison & Steele.  */
  union
  {
    long l;
    char c[sizeof (long)];
  } u;
  u.l = 1;
  return (u.c[sizeof (long) - 1] == 1);
}

