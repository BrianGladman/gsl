#include <config.h>
#include <stdio.h>
#include <math.h>
#include <gsl_ieee_utils.h>

#include "ieee_utils.h"

/* A table of sign characters, 0=positive, 1=negative. We print a space
   instead of a unary + sign for compatibility with bc */

static char signs[2]={' ','-'} ;  

/* A table of character representations of nybbles */

static char nybble[16][4]={
  "0000", "0001", "0010", "0011",
  "0100", "0101", "0110", "0111",
  "1000", "1001", "1010", "1011",
  "1100", "1101", "1110", "1111"
}  ;

void 
gsl_ieee_printf_float (const float * x) {
  gsl_ieee_float_rep r ;
  gsl_ieee_float_to_rep(x, &r) ;

  /* output is compatible with bc (with ibase=2), mant*2^expb */  

  switch (r.type)
    {
    case GSL_IEEE_TYPE_NAN:
      printf("NaN") ;
      break ;
    case GSL_IEEE_TYPE_INF:
      printf("%cInf", signs[r.sign]) ;
      break ;
    case GSL_IEEE_TYPE_NORMAL:
      printf("%c1.%.23s*2^%d", signs[r.sign], r.bits, r.exponent) ;
      break ;
    case GSL_IEEE_TYPE_DENORMAL:
      printf("%c0.%.23s*2^%d", signs[r.sign], r.bits, r.exponent + 1) ;
      break ;
    case GSL_IEEE_TYPE_ZERO:
      printf("%c0", signs[r.sign]) ;
      break ;
    default:
      printf("[non-standard IEEE float]") ;
    }
}

void
gsl_ieee_printf_double (const double * x) {
  gsl_ieee_double_rep r ;
  gsl_ieee_double_to_rep (x, &r) ;

  /* output is compatible with bc (with ibase=2), mant*2^expb */  

  switch (r.type)
    {
    case GSL_IEEE_TYPE_NAN:
      printf("NaN") ;
      break ;
    case GSL_IEEE_TYPE_INF:
      printf("%cInf", signs[r.sign]) ;
      break ;
    case GSL_IEEE_TYPE_NORMAL:
      printf("%c1.%.52s*2^%d", signs[r.sign], r.bits, r.exponent) ;
      break ;
    case GSL_IEEE_TYPE_DENORMAL:
      printf("%c0.%.52s*2^%d", signs[r.sign], r.bits, r.exponent + 1) ;
      break ;
    case GSL_IEEE_TYPE_ZERO:
      printf("%c0", signs[r.sign]) ;
      break ;
    default:
      printf("[non-standard IEEE double]") ;
    }
}

/* For the IEEE float format the bits are found from the following
   masks,
   
   sign      = 0x80000000  
   exponent  = 0x7f800000 
   mantisssa = 0x007fffff  

   For the IEEE double format the masks are,

   sign      = 0x8000000000000000  
   exponent  = 0x7ff0000000000000 
   mantissa  = 0x000fffffffffffff
   
   */

void 
gsl_ieee_float_to_rep (const float * x, gsl_ieee_float_rep * r)
{
  int e, non_zero;

  union { 
    float f;
    struct  { 
      unsigned char byte[4] ;
    } ieee ;
  } u;
  
  u.f = *x ; 

  if (little_endian_p())
    make_float_bigendian(&(u.f)) ;
  
  r->sign = u.ieee.byte[3]>>7 ;

  e = (u.ieee.byte[3] & 0x7f) << 1 | (u.ieee.byte[2] & 0x80)>>7 ; 
  
  r->exponent = e - 127 ;

  sprint_byte((u.ieee.byte[2] & 0x7f) << 1,r->bits) ;
  sprint_byte(u.ieee.byte[1],r->bits + 7) ;
  sprint_byte(u.ieee.byte[0],r->bits + 15) ;

  non_zero = u.ieee.byte[0] || u.ieee.byte[1] || (u.ieee.byte[2] & 0x7f);

  r->type = determine_ieee_type (e, non_zero, 255) ;
}

void 
gsl_ieee_double_to_rep (const double * x, gsl_ieee_double_rep * r)
{

  int e, non_zero;

  union 
  { 
    double d;
    struct  { 
      unsigned char byte[8];
    } ieee ;
  } u;

  u.d= *x ; 
  
  if (little_endian_p())
    make_double_bigendian(&(u.d)) ;
  
  r->sign = u.ieee.byte[7]>>7 ;

  e =(u.ieee.byte[7] & 0x7f)<<4 ^ (u.ieee.byte[6] & 0xf0)>>4 ;
  
  r->exponent = e - 1023 ;

  sprint_nybble(u.ieee.byte[6],r->bits) ;
  sprint_byte(u.ieee.byte[5],r->bits + 4) ;
  sprint_byte(u.ieee.byte[4],r->bits + 12) ;
  sprint_byte(u.ieee.byte[3],r->bits + 20) ; 
  sprint_byte(u.ieee.byte[2],r->bits + 28) ;
  sprint_byte(u.ieee.byte[1],r->bits + 36) ;
  sprint_byte(u.ieee.byte[0],r->bits + 44) ;

  non_zero = (u.ieee.byte[0] || u.ieee.byte[1] || u.ieee.byte[2]
	      || u.ieee.byte[3] || u.ieee.byte[4] || u.ieee.byte[5] 
	      || (u.ieee.byte[6] & 0x0f)) ;

  r->type = determine_ieee_type (e, non_zero, 2047) ;
}
	  
void
sprint_nybble(int i, char *s)
{
  char *c ;
  c=nybble[i & 0x0f ];
  *s=c[0] ;  *(s+1)=c[1] ;  *(s+2)=c[2] ;  *(s+3)=c[3] ;
} 

void
sprint_byte(int i, char *s)
{
  char *c ;
  c=nybble[(i & 0xf0)>>4];
  *s=c[0] ;  *(s+1)=c[1] ;  *(s+2)=c[2] ;  *(s+3)=c[3] ;
  c=nybble[i & 0x0f];
  *(s+4)=c[0] ;  *(s+5)=c[1] ;  *(s+6)=c[2] ;  *(s+7)=c[3] ;
} 

int 
determine_ieee_type (int exponent, int max_exponent, int non_zero)
{
  if (exponent == max_exponent)
    {
      if (non_zero)
	{
	  return GSL_IEEE_TYPE_NAN ;
	}
      else
	{
	  return GSL_IEEE_TYPE_INF ;
	}
    }
  else if (exponent == 0)
    {
      if (non_zero)
	{
	  return GSL_IEEE_TYPE_DENORMAL ;
	}
      else
	{
	  return GSL_IEEE_TYPE_ZERO ;
	}
    }
  else
    {
      return GSL_IEEE_TYPE_NORMAL ;
    }
}



