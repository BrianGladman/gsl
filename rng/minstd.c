#include <config.h>
#include <stdlib.h>
#include <gsl_rng.h>

/* MINSTD is Park and Miller's minimal standard generator (i.e. it's
   not particularly good).

   The sequence is

   x_{n+1} = (a x_n) mod m

   with a = 16807 and m = 2^31 - 1 = 2147483647. The seed specifies
   the initial value, x_1.  

   The theoretical value of x_{10001} is 1043618065, starting with a
   seed of x_1 = 1.

   The period of this generator is 2^31.

   From: Park and Miller, "Random Number Generators: Good ones are
   hard to find" Communications of the ACM, October 1988, Volume 31,
   No 10, pages 1192-1201. */

unsigned long int minstd_get (void *vstate);
void minstd_set (void *state, unsigned long int s);

static const long int m = 2147483647, a = 16807, q = 127773, r = 2836;

typedef struct
  {
    unsigned long int x;
  }
minstd_state_t;

unsigned long int
minstd_get (void *vstate)
{
  minstd_state_t *state = (minstd_state_t *) vstate;

  const unsigned long int x = state->x;

  const long int h = x / q;
  const long int t = a * (x - h * q) - h * r;

  if (t < 0)
    {
      state->x = t + m;
    }
  else
    {
      state->x = t;
    }

  return state->x;
}

void
minstd_set (void *vstate, unsigned long int s)
{
  minstd_state_t *state = (minstd_state_t *) vstate;

  if (s == 0)
    s = 1;	/* default seed is 1 */

  state->x = s;

  return;
}

static const gsl_rng_type minstd_type =
{"minstd",			/* name */
 2147483646,			/* RAND_MAX */
 0,          			/* RAND_MIN */
 sizeof (minstd_state_t),
 &minstd_set,
 &minstd_get};

const gsl_rng_type *gsl_rng_minstd = &minstd_type;
