#include <config.h>
#include <stdlib.h>
#include <gsl_ieee_utils.h>

#include "ieee_utils.h"

int 
little_endian_p (void) {
  /* Are we little or big endian?  From Harbison & Steele.  */
  union
  {
    long l;
    char c[sizeof (long)];
  } u;
  u.l = 1;
  return (u.c[sizeof (long) - 1] == 1);
}

int 
endianness (void) {
  /* Determine true endianness, big endian 4321, little endian 1234. */
  union
  {
    long l;
    char c[sizeof (long)];
  } u;
  u.l = 0x04030201;
  return (u.c[0]+10*(u.c[1]+10*(u.c[2]+10*u.c[3])));
}

int 
double_endianness (void) {
  /* Determine true of doubles endianness, big endian 4321, little
     endian 1234. */
  int i, w ;
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
  
  return (w);
}

void setup_dynamic_endianness(int *b0,int *b1,int *b2,int *b3,
			      int *b4,int *b5,int *b6,int *b7)
{

  int i, j, flag=0;

  union
    {
      double d ;
      char c[sizeof(double)];
    } u;

  /* this corresponds to 0x0706050403020100 in big endian format */

  u.d=7.9499288951273625e-275 ;

  *b0=u.c[0] ;  *b1=u.c[1] ;  *b2=u.c[2] ;  *b3=u.c[3] ;
  *b4=u.c[4] ;  *b5=u.c[5] ;  *b6=u.c[6] ;  *b7=u.c[7] ;


  /* check for any out of range values , or corruption of the 00 01 02
     03 04 05 06 07 permutation */

  for(i=0 ; i<8 ; i++)
    {
      if(u.c[i]>7) flag=1 ;
    }

  for(i=0 ; i<7 ; i++)
    {
      for(j=i+1 ; j<8 ; j++)
	{
	  if(u.c[i]==u.c[j]) flag=1 ;
	}
    }
  
  if(flag)
    {
      exit(EXIT_FAILURE) ;
    } ;
}
