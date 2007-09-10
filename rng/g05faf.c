/* rng/g05faf.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 James Theiler, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <config.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>

/* This is the NAG G05FAF generator. The generator returns the
   upper 32 bits from each term of the sequence,

   x_{n+1} = (a x_n) mod m 

   using 59-bit unsigned arithmetic, with a = 13^{13} and m =
   2^59. The seed specifies the upper 32 bits of the initial value,
   x_1, with the lower 16 bits set to 0x330E.

   The theoretical value of x_{10001} is 244131582646046.

   The period of this generator is ? FIXME (probably around 2^48). */

static inline void g05faf_advance (void *vstate);
static unsigned long int g05faf_get (void *vstate);
static double g05faf_get_double (void *vstate);
static void g05faf_set (void *state, unsigned long int s);

static const unsigned short int a0 = 0xFD ;
static const unsigned short int a1 = 0xC5 ;
static const unsigned short int a2 = 0x23 ;
static const unsigned short int a3 = 0x9B ;
static const unsigned short int a4 = 0x76 ;
static const unsigned short int a5 = 0x13 ;
static const unsigned short int a6 = 0x01 ;
static const unsigned short int a7 = 0x00 ;

typedef struct
  {
    unsigned short int x0, x1, x2, x3, x4, x5, x6, x7 ;
  }
g05faf_state_t;

static inline void
g05faf_advance (void *vstate)
{
  g05faf_state_t *state = (g05faf_state_t *) vstate;

  const unsigned short int x0 = state->x0 ;
  const unsigned short int x1 = state->x1 ;
  const unsigned short int x2 = state->x2 ;
  const unsigned short int x3 = state->x3 ;
  const unsigned short int x4 = state->x4 ;
  const unsigned short int x5 = state->x5 ;
  const unsigned short int x6 = state->x6 ;
  const unsigned short int x7 = state->x7 ;

  unsigned long int a ;

  /* This looks like it will be pretty slow. Maybe someone can produce
     a more efficient implementation (it has to be portable though) */
  
  a = a0 * x0 ;
  state->x0 = (a & 0x000000FFUL) ;
 
  a >>= 8 ;
  a += a0 * x1 + a1 * x0 ;
  state->x1 = (a & 0x000000FFUL) ;
  
  a >>= 8 ;
  a += a0 * x2 + a1 * x1 + a2 * x0 ;
  state->x2 = (a & 0x000000FFUL) ;

  a >>= 8 ;
  a += a0 * x3 + a1 * x2 + a2 * x1 + a3 * x0 ;
  state->x3 = (a & 0x000000FFUL) ;

  a >>= 8 ;
  a += a0 * x4 + a1 * x3 + a2 * x2 + a3 * x1 + a4 * x0 ;
  state->x4 = (a & 0x000000FFUL) ;

  a >>= 8 ;
  a += a0 * x5 + a1 * x4 + a2 * x3 + a3 * x2 + a4 * x1 + a5 * x0 ;
  state->x5 = (a & 0x000000FFUL) ;

  a >>= 8 ;
  a += a0 * x6 + a1 * x5 + a2 * x4 + a3 * x3 + a4 * x2 + a5 * x1 + a6 * x0 ;
  state->x6 = (a & 0x000000FFUL) ;

  a >>= 8 ;
  a += (a0 * x7 + a1 * x6 + a2 * x5 + a3 * x4 + a4 * x3 + a5 * x2 + a6 * x1 
        + a7 * x0) ;
  state->x7 = (a & 0x00000007UL) ;

}

static unsigned long int 
g05faf_get (void *vstate)
{
  g05faf_state_t *state = (g05faf_state_t *) vstate;

  g05faf_advance (state) ;

  return ((state->x7 << 29) + (state->x6 << 21) +  
          (state->x5 << 13) + (state->x4 << 5) + (state->x3 >> 3));
}

static double
g05faf_get_double (void * vstate)
{
  g05faf_state_t *state = (g05faf_state_t *) vstate;

  g05faf_advance (state) ;  

  return (ldexp((double) state->x7, -3)
          + ldexp((double) state->x6, -11) 
          + ldexp((double) state->x5, -19)
          + ldexp((double) state->x4, -27)
          + ldexp((double) state->x3, -35)
          + ldexp((double) state->x2, -43)
          + ldexp((double) state->x1, -51)
          + ldexp((double) state->x0, -59)) ;
}

static void
g05faf_set (void *vstate, unsigned long int s)
{
  g05faf_state_t *state = (g05faf_state_t *) vstate;

  if (s == 0)  /* default seed */
    s = 1 ;
  
  /* FIXME: I have no idea what the nag seeding procedure is */

  return;
}

static const gsl_rng_type g05faf_type =
{"g05faf",                      /* name */
 0xffffffffUL,                  /* RAND_MAX */
 0,                             /* RAND_MIN */
 sizeof (g05faf_state_t),
 &g05faf_set,
 &g05faf_get,
 &g05faf_get_double
};

const gsl_rng_type *gsl_rng_g05faf = &g05faf_type;
