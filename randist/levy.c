/* randist/levy.c
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
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

/* The stable Levy probability distributions have the form

   p(x) dx = (1/(2 pi)) \int dt exp(it(mu-x) - |ct|^alpha)

   with 0 < alpha <= 2

   For alpha = 1, we get the Cauchy distribution
   For alpha = 2, we get the Gaussian distribution

   I found this in Chapter 5 of Bratley, Fox and Schrage "A Guide to
   Simulation". The original reference was,

   J.M. Chambers, C.L. Mallows and B. W. Stuck. "A method for
   simulating stable random variates". Journal of the American
   Statistical Association, JASA 71 340-344 (1976).

   */

double
gsl_ran_levy (const gsl_rng * r, const double mu, const double a)
{
  double u, v, t, s;

  do
    {
      u = M_PI * (gsl_rng_uniform_pos (r) - 0.5);
    }
  while (u == 0);
  
  do
    {
      v = gsl_ran_exponential (r, 1.0);
    } 
  while (v == 0);
  
  t = sin(a * u) / pow(cos(u), 1/a) ;

  if (a == 1)
    {
      s = 1 ;
    }
  else 
    {
      s = pow(cos((1 - a)*u) / v, (1-a)/a) ;
    }
  
  return mu * t * s;
}

