/* rng/default.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 James Theiler, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include <config.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_errno.h>

/* The initial defaults are defined in the file mt.c, so we can get
   access to the static parts of the default generator. */

static void
  check (const gsl_rng_type ** def, const gsl_rng_type * test, const char *p);

const gsl_rng_type *
gsl_rng_env_setup (void)
{
  unsigned long int seed = 0;
  const char *p = getenv ("GSL_RNG_TYPE");

  if (p)
    {
      gsl_rng_default = 0;

      /* check GSL_RNG_TYPE against the names of all the generators */

      check (&gsl_rng_default, gsl_rng_slatec, p);
      check (&gsl_rng_default, gsl_rng_cmrg, p);
      check (&gsl_rng_default, gsl_rng_gfsr4, p);
      check (&gsl_rng_default, gsl_rng_minstd, p);
      check (&gsl_rng_default, gsl_rng_mrg, p);
      check (&gsl_rng_default, gsl_rng_mt19937, p);
      check (&gsl_rng_default, gsl_rng_r250, p);
      check (&gsl_rng_default, gsl_rng_ran0, p);
      check (&gsl_rng_default, gsl_rng_ran1, p);
      check (&gsl_rng_default, gsl_rng_ran2, p);
      check (&gsl_rng_default, gsl_rng_ran3, p);
      check (&gsl_rng_default, gsl_rng_rand, p);
      check (&gsl_rng_default, gsl_rng_rand48, p);
      check (&gsl_rng_default, gsl_rng_random8_bsd, p);
      check (&gsl_rng_default, gsl_rng_random8_glibc2, p);
      check (&gsl_rng_default, gsl_rng_random8_libc5, p);
      check (&gsl_rng_default, gsl_rng_random128_bsd, p);
      check (&gsl_rng_default, gsl_rng_random128_glibc2, p);
      check (&gsl_rng_default, gsl_rng_random128_libc5, p);
      check (&gsl_rng_default, gsl_rng_random256_bsd, p);
      check (&gsl_rng_default, gsl_rng_random256_glibc2, p);
      check (&gsl_rng_default, gsl_rng_random256_libc5, p);
      check (&gsl_rng_default, gsl_rng_random32_bsd, p);
      check (&gsl_rng_default, gsl_rng_random32_glibc2, p);
      check (&gsl_rng_default, gsl_rng_random32_libc5, p);
      check (&gsl_rng_default, gsl_rng_random64_bsd, p);
      check (&gsl_rng_default, gsl_rng_random64_glibc2, p);
      check (&gsl_rng_default, gsl_rng_random64_libc5, p);
      check (&gsl_rng_default, gsl_rng_random_bsd, p);
      check (&gsl_rng_default, gsl_rng_random_glibc2, p);
      check (&gsl_rng_default, gsl_rng_random_libc5, p);
      check (&gsl_rng_default, gsl_rng_randu, p);
      check (&gsl_rng_default, gsl_rng_ranf, p);
      check (&gsl_rng_default, gsl_rng_ranlux, p);
      check (&gsl_rng_default, gsl_rng_ranlux389, p);
      check (&gsl_rng_default, gsl_rng_ranlxd1, p);
      check (&gsl_rng_default, gsl_rng_ranlxd2, p);
      check (&gsl_rng_default, gsl_rng_ranmar, p);
      check (&gsl_rng_default, gsl_rng_taus, p);
      check (&gsl_rng_default, gsl_rng_transputer, p);
      check (&gsl_rng_default, gsl_rng_tt800, p);
      check (&gsl_rng_default, gsl_rng_uni, p);
      check (&gsl_rng_default, gsl_rng_uni32, p);
      check (&gsl_rng_default, gsl_rng_vax, p);
      check (&gsl_rng_default, gsl_rng_zuf, p);

      if (gsl_rng_default == 0)
	{
	  GSL_ERROR_RETURN ("unknown generator", GSL_EINVAL, 0);
	}

      fprintf (stderr, "GSL_RNG_TYPE=%s\n", gsl_rng_default->name);
    }
  else
    {
      gsl_rng_default = gsl_rng_mt19937;
    }

  p = getenv ("GSL_RNG_SEED");

  if (p)
    {
      seed = strtoul (p, 0, 0);
      fprintf (stderr, "GSL_RNG_SEED=%lu\n", seed);
    };

  gsl_rng_default_seed = seed;

  return gsl_rng_default;
}


static void
check (const gsl_rng_type ** def, const gsl_rng_type * test, const char *p)
{
  if (*def)
    return;

  if (strcmp (p, test->name) == 0)
    {
      *def = test;
    }
}
