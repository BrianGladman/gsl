#include <config.h>
#include <stdlib.h>
#include <gsl_rng.h>

/* This is the INMOS Transputer Development System generator. The sequence is,

   x_{n+1} = (a x_n) mod m

   with a = 1664525 and m = 2^32. The seed specifies the initial
   value, x_1.

   The theoretical value of x_{10001} is FIXME.

   The period of this generator is 2^32. */

unsigned long int tds_get (void *vstate);
void tds_set (void *state, unsigned long int s);

typedef struct
  {
    unsigned long int x;
  }
tds_state_t;

unsigned long int
tds_get (void *vstate)
{
  tds_state_t *state = (tds_state_t *) vstate;

  state->x = (1664525 * state->x) & 0xffffffffUL;

  return state->x;
}

void
tds_set (void *vstate, unsigned long int s)
{
  tds_state_t *state = (tds_state_t *) vstate;

  if (s == 0)
    s = 1 ;   /* default seed is 1. */

  state->x = s;

  return;
}

static const gsl_rng_type tds_type =
{"tds",				/* name */
 0xffffffffUL,			/* RAND_MAX */
 1,				/* RAND_MIN */
 sizeof (tds_state_t),
 &tds_set,
 &tds_get};

const gsl_rng_type *gsl_rng_tds = &tds_type;
