/* This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License as
   published by the Free Software Foundation; either version 2 of the
   License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.  You should have received
   a copy of the GNU General Public License along with this program;
   if not, write to the Free Foundation, Inc., 59 Temple Place, Suite
   330, Boston, MA 02111-1307 USA

   Copyright (C) 1998 Brian Gough, reimplemented for GSL

   Original work was copyright (C) 1997 Makoto Matsumoto and Takuji
   Nishimura. Coded by Takuji Nishimura, considering the suggestions
   by Topher Cooper and Marc Rieffel in July-Aug. 1997, "A C-program
   for MT19937: Integer version (1998/4/6)"

   When you use this, send an email to: matumoto@math.keio.ac.jp with
   an appropriate reference to your work.  */

#include <gsl_rng.h>

unsigned long int mt_get (void *vstate);
void mt_set (void *state, unsigned long int s);

#define N 624	/* Period parameters */
#define M 397
const unsigned long UPPER_MASK = 0x80000000UL;	/* most significant w-r bits */
const unsigned long LOWER_MASK = 0x7fffffffUL;	/* least significant r bits */

typedef struct
  {
    int mti;
    unsigned long mt[N];
  }
mt_state_t;

unsigned long
mt_get (void *vstate)
{
  mt_state_t *state = (mt_state_t *) vstate;

  unsigned long y;
  const unsigned long int mag01[2] =
  {0x00000000UL, 0x9908b0dfUL};
  unsigned long int *const mt = state->mt;

  if (state->mti >= N)
    {	/* generate N words at one time */
      int kk;

      for (kk = 0; kk < N - M; kk++)
	{
	  y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
	  mt[kk] = mt[kk + M] ^ (y >> 1) ^ mag01[y & 0x1];
	}
      for (; kk < N - 1; kk++)
	{
	  y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
	  mt[kk] = mt[kk + (M - N)] ^ (y >> 1) ^ mag01[y & 0x1];
	}
      y = (mt[N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
      mt[N - 1] = mt[M - 1] ^ (y >> 1) ^ mag01[y & 0x1];

      state->mti = 0;
    }

  /* Tempering */

  y = mt[state->mti];
  y ^= (y >> 11);
  y ^= (y << 7) & 0x9d2c5680UL;
  y ^= (y << 15) & 0xefc60000UL;
  y ^= (y >> 18);

  state->mti++;

  return y;
}

void
mt_set (void *vstate, unsigned long int s)
{
  mt_state_t *state = (mt_state_t *) vstate;
  int i;

  if (s == 0)
    s = 4357;	/* the default seed is 4357 */

  state->mt[0] = s & 0xffffffffUL;

#define LCG(n) ((69069 * n) & 0xffffffffUL)

  for (i = 1; i < N; i++)
    state->mt[i] = LCG (state->mt[i - 1]);

  state->mti = i;
}

static const gsl_rng_type mt_type =
{"mt19937",			/* name */
 0xffffffffUL,			/* RAND_MAX  */
 sizeof (mt_state_t),
 &mt_set,
 &mt_get};

const gsl_rng_type *gsl_rng_mt19937 = &mt_type;
