/* sum/demo.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
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
#include <gsl/gsl_math.h>
#include <gsl/gsl_sum.h>

#define N 20

int 
main (void)
{
  double t[N];
  double sum_accel, err;
  double sum = 0;
  int n;
  
  gsl_sum_levin_u_workspace * w = gsl_sum_levin_u_alloc (N);

  const double zeta_2 = M_PI * M_PI / 6.0;
  
  /* terms for zeta(2) = \sum_{n=1}^{\infty} 1/n^2 */

  for (n = 0; n < N; n++)
    {
      double np1 = n + 1.0;
      t[n] = 1.0 / (np1 * np1);
      sum += t[n] ;
    }
  
  gsl_sum_levin_u_accel (t, N, w, &sum_accel, &err);

  printf("term-by-term sum = % .16f using %d terms\n", sum, N) ;

  printf("term-by-term sum = % .16f using %d terms\n", 
         w->sum_plain, w->terms_used) ;

  printf("exact value      = % .16f\n", zeta_2) ;
  printf("accelerated sum  = % .16f using %d terms\n", 
         sum_accel, w->terms_used) ;

  printf("estimated error  = % .16f\n", err) ;
  printf("actual error     = % .16f\n", sum_accel - zeta_2) ;

  gsl_sum_levin_u_free (w);

  return 0;
}
