#include <config.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>

/* This is a reincarnation of the infamously bad RANDU generator.
   The sequence is,

   x_{n+1} = (a x_n) mod m

   with a = 65539 and m = 2^31 = 2147483648. The seed specifies
   the initial value, x_1.

   The theoretical value of x_{10001} is 1623524161.

   The period of this generator is 2^29.

   Note: Knuth describes this generator as "really horrible". 

   From: Park and Miller, "Random Number Generators: Good ones are
   hard to find" Communications of the ACM, October 1988, Volume 31,
   No 10, pages 1192-1201. */

static inline unsigned long int randu_get (void *vstate);
static double randu_get_double (void *vstate);
static void randu_set (void *state, unsigned long int s);

static const long int a = 65539;
static const unsigned long int m = 2147483648UL;

typedef struct
  {
    unsigned long int x;
  }
randu_state_t;

static inline unsigned long int
randu_get (void *vstate)
{
  randu_state_t *state = (randu_state_t *) vstate;

  /* The following line relies on unsigned 32-bit arithmetic */

  state->x = (a * state->x) & 0x7fffffffUL;

  return state->x;
}

static double
randu_get_double (void *vstate)
{
  return randu_get (vstate) / 2147483648.0 ;
}

static void
randu_set (void *vstate, unsigned long int s)
{
  randu_state_t *state = (randu_state_t *) vstate;

  if (s == 0)
    s = 1;	/* default seed is 1 */

  state->x = s;

  return;
}

static const gsl_rng_type randu_type =
{"randu",			/* name */
 0x7fffffffUL,			/* RAND_MAX */
 1,       			/* RAND_MIN */
 sizeof (randu_state_t),
 &randu_set,
 &randu_get,
 &randu_get_double};

const gsl_rng_type *gsl_rng_randu = &randu_type;
