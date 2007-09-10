/* ieee-utils/demo.c
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

#include <float.h>
#include <stdio.h>
#include <gsl/gsl_ieee_utils.h>

int
main () 
{
  int i ;
  double x = 10 ;
  float f = 1.0/3.0 ;
  double d = 1.0/3.0 ;

  double fd = f ; /* promote from float to double */

  gsl_ieee_env_setup() ;
  
  printf("     float 1/3 = ") ; gsl_ieee_printf_float(&f) ; printf("\n") ;
  printf("promoted float = ") ; gsl_ieee_printf_double(&fd) ; printf("\n") ;
  printf("    double 1/3 = ") ; gsl_ieee_printf_double(&d) ; printf("\n") ;
  
  for (i=0;i<10;i++) {
    x = 0.5 *(x + 2.0/x) ;
    printf("%.18g ",x) ;
    gsl_ieee_printf_double(&x) ; 
    printf("\n") ;
  }

  x = 10 * x * GSL_DBL_MAX ;

  printf("%.18g ",x) ; gsl_ieee_printf_double(&x) ; printf("\n") ;

  x = x / 0 ; ;

  printf("%.18g ",x) ; gsl_ieee_printf_double(&x) ; printf("\n") ;
  
  f = -1.0/3.0 ;

  while (f < 0) {
    f = f / 10 ;
    printf("%.18g ",f) ; gsl_ieee_printf_float(&f) ; printf("\n") ;
  }

}

