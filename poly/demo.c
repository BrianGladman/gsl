/* poly/demo.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
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
#include <gsl/gsl_poly.h>

int
main ()
{
  int i;
  double a[6] = { -1, 0, 0, 0, 0, 1 };  /*  P(x) =  x^5 - 1  */
  double z[10];

  gsl_poly_complex_workspace * w = gsl_poly_complex_workspace_alloc (6) ;
  
  gsl_poly_complex_solve (a, 6, w, z) ;

  gsl_poly_complex_workspace_free (w) ;

  for (i = 0; i < 5 ; i++)
    {
      printf("z%d = %+.18f %+.18f\n", i, z[2*i], z[2*i+1]) ;
    }
}
