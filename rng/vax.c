#include <config.h>
#include <stdlib.h>
#include <gsl_rng.h>

/* This is the old vax generator MTH$RANDOM. The sequence is,

   x_{n+1} = (a x_n + c) mod m

   with a = 69069, c = 1 and m = 2^32. The seed specifies the initial
   value, x_1.

   The theoretical value of x_{10001} is 3051034865.

   The period of this generator is 2^32. */

unsigned long int vax_get (void *vstate);
void vax_set (void *state, unsigned long int s);

typedef struct
  {
    unsigned long int x;
  }
vax_state_t;

unsigned long int
vax_get (void *vstate)
{
  vax_state_t *state = (vax_state_t *) vstate;

  state->x = (69069 * state->x + 1) & 0xffffffffUL;

  return state->x;
}

void
vax_set (void *vstate, unsigned long int s)
{
  vax_state_t *state = (vax_state_t *) vstate;

  /* default seed is 0. The constant term c stops the series from
     collapsing to 0,0,0,0,0,... */

  state->x = s;

  return;
}

static const gsl_rng_type vax_type =
{"vax",				/* name */
 0xffffffffUL,			/* RAND_MAX */
 0,				/* RAND_MIN */
 sizeof (vax_state_t),
 &vax_set,
 &vax_get};

const gsl_rng_type *gsl_rng_vax = &vax_type;
