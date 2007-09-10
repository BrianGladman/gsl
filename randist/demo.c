/* randist/demo.c
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

#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

int
main ()
{
  gsl_rng * r ;

  int i, n = 10;
  double mu = 3.0;

  /* create a generator chosen by the environment variable GSL_RNG_TYPE */

  gsl_rng_env_setup();
  
  r = gsl_rng_alloc (gsl_rng_default);

  /* print n random variates chosen from the poisson distribution with
     mean parameter mu */

  for (i = 0; i < n; i++) 
    {
      unsigned int k = gsl_ran_poisson (r, mu);
      
      printf(" %u", k);
    }

  printf("\n");
}
