#include <stdio.h>
#include <gsl_ieee_utils.h>

int 
main (void)
{
  float x = 1.000 ;
  double y = -8.000 ;
  float d = 0.5 ;
  int i, w ;

  y = 2.1 ; printf("%g ",y) ;    gsl_ieee_printf_double(&y) ; printf("\n") ;

  for (i = 1; i<30; i++) {
    printf("%g ",x) ;    gsl_ieee_printf_float(&x) ; printf("\n") ;
    x += d ;
    d = d / 2 ;
  }

  y = 1 ;
  d = 0.5 ;
 
  for (i = 1; i<400; i++) {
    printf("%g ",y) ;    gsl_ieee_printf_double(&y) ; printf("\n") ;
    y *= 2 ;
  }


  {
    union
    {
      double d ;
      char c[sizeof(double)];
    } u;
    
    /* this corresponds to 0x0706050403020100 in big endian format */
    u.d=7.9499288951273625e-275 ;
    
    w=u.c[7]+1 ;
    for(i=6 ; i>=0 ; i--)
      {
	w=10*w+u.c[i]+1 ;
      } ;
    
    printf("%d\n",w);
  }

  {
    union
    {
      long l;
      char c[sizeof (long)];
    } u;
    u.l = 0x04030201;
    printf("%d\n",u.c[0]+10*(u.c[1]+10*(u.c[2]+10*u.c[3]))) ;
  }
  return 0 ;
}

