#include <float.h>
#include <stdio.h>
#include <gsl_ieee_utils.h>

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

  x = 10 * x * DBL_MAX ;

  printf("%.18g ",x) ; gsl_ieee_printf_double(&x) ; printf("\n") ;

  x = x / 0 ; ;

  printf("%.18g ",x) ; gsl_ieee_printf_double(&x) ; printf("\n") ;
}

