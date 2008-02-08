/* Modified from the following source for GSL by James Scott */ 
/* Globals have been removed, but now everything is in function
   parameters... could do with tidying up */ 

/*===================================================================*/
/* C program for distribution from the Combinatorial Object Server.  */
/* Program to enumerate all irreducible and/or primitive polynomials */
/* over GF(2).  This is the same version described in the book       */
/* "Combinatorial Generation", Frank Ruskey, to appear.              */
/* The program can be modified, translated to other languages, etc., */
/* so long as proper acknowledgement is given (author and source).   */
/* Not to be used for commercial purposes without prior consent.     */
/* The latest version of this program may be found at the site       */
/* http://www.theory.csc.uvic.ca/~cos/inf/neck/Polynomial.html       */
/* Copyright 1997,1998 F. Ruskey and K. Cattell                      */
/*===================================================================*/
                                                                            
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>



/* 
   typedef unsigned long long coslong; 
   #define MAX 64
*/ 


/* For GSL we're only interested in polynomials representable in an 'int' */ 
/* MAX should be no. bits in 'coslong' */ 
typedef  int coslong; 
#define MAX 32 


#define ONE ((coslong)1)
#define TWO ((coslong)2)



struct pair {
  coslong poly;      /* some primitive polynomial */
  coslong pow2m1;    /* 2^n - 1 */
};

static const struct pair poly_table[MAX+1] = {  /* prime factorization of 2^n-1 */
  /*  0 */  {  1,          0U },  /* (0)                               */
  /*  1 */  {  1,          1U },  /* (1)                               */
  /*  2 */  {  3,          3U },  /* (3)          it's prime!          */
  /*  3 */  {  3,          7U },  /* (7)          it's prime!          */
  /*  4 */  {  3,         15U },  /* (3) (5)                           */
  /*  5 */  {  5,         31U },  /* (31)         it's prime!          */
  /*  6 */  {  3,         63U },  /* (3) (3) (7)                       */
  /*  7 */  {  3,        127U },  /* (127)        it's prime!          */
  /*  8 */  { 29,        255U },  /* (3) (5) (17)                      */
  /*  9 */  { 17,        511U },  /* (7) (73)                          */
  /* 10 */  {  9,       1023U },  /* (3) (11) (31)                     */
  /* 11 */  {  5,       2047U },  /* (23) (89)                         */
  /* 12 */  { 83,       4095U },  /* (3) (3) (5) (7) (13)              */
  /* 13 */  { 27,       8191U },  /* (8191)       it's prime!          */
  /* 14 */  { 43,      16383U },  /* (3) (43) (127)                    */
  /* 15 */  {  3,      32767U },  /* (7) (31) (151)                    */
  /* 16 */  { 45,      65535U },  /* (3) (5) (17) (257)                */
  /* 17 */  {  9,     131071U },  /* (131071)     it's prime!          */
  /* 18 */  { 39,     262143U },  /* (3) (3) (7) (19) (73)             */
  /* 19 */  { 39,     524287U },  /* (524287)     it's prime!          */
  /* 20 */  {  9,    1048575U },  /* (3) (5) (5) (11) (31) (41)        */
  /* 21 */  {  5,    2097151U },  /* (7) (7) (127) (337)               */
  /* 22 */  {  3,    4194303U },  /* (3) (23) (89) (683)               */
  /* 23 */  { 33,    8388607U },  /* (47) (178481)                     */
  /* 24 */  { 27,   16777215U },  /* (3) (3) (5) (7) (13) (17) (241)   */
  /* 25 */  {  9,   33554431U },  /* (31) (601) (1801)                 */
  /* 26 */  { 71,   67108863U },  /* (3) (8191) (2731)                 */
  /* 27 */  { 39,  134217727U },  /* (7) (73) (262657)                 */
  /* 28 */  {  9,  268435455U },  /* (3) (5) (29) (43) (113) (127)     */
  /* 29 */  {  5,  536870911U },  /* (233) (1103) (2089)               */
  /* 30 */  { 83, 1073741823U },  /* (3) (3) (7) (11) (31) (151) (331) */
  /* 31 */  {  9, 2147483647U },  /* (2147483647) it's prime!          */
  /* 32 */  {175, 4294967295U }   /* (3) (5) (17) (257) (65537)        */
#if (MAX > 32) 
  ,
  /* 33 */  { 83, 8589934591U },  /* (7) (23) (89) (599479)            */
  /* 34 */  {231, 17179869183U        },/*131071.3.43691*/
  /* 35 */  {5,   34359738367U        },/*71.122921.31.127*/
  /* 36 */  {119, 68719476735U        },/*7.73.3.3.19.13.5.37.109*/
  /* 37 */  {63,  137438953471U       },/*616318177.223*/
  /* 38 */  {99,  274877906943U       },/*524287.3.174763*/
  /* 39 */  {17,  549755813887U       },/*7.8191.121369.79*/
  /* 40 */  {57,  1099511627775U      },/*31.3.11.5.5.41.17.61681*/
  /* 41 */  {9,   2199023255551U      },/*164511353.13367*/
  /* 42 */  {63,  4398046511103U      },/*337.7.7.127.43.3.3.5419*/
  /* 43 */  {89,  8796093022207U      },/*2099863.431.9719*/
  /* 44 */  {101, 17592186044415U     },/*23.89.3.683.5.397.2113*/
  /* 45 */  {27,  35184372088831U     },/*7.151.73.31.631.23311*/
  /* 46 */  {303, 70368744177663U     },/*178481.47.3.2796203*/
  /* 47 */  {33,  140737488355327U    },/*2351.13264529.4513*/
  /* 48 */  {183, 281474976710655U    },/*7.3.3.13.5.17.241.257.97.673*/
  /* 49 */  {113, 562949953421311U    },/*127.4432676798593*/
  /* 50 */  {29,  1125899906842623U   },/*601.1801.31.3.11.251.4051*/
  /* 51 */  {75,  2251799813685247U   },/*7.131071.11119.2143.103*/
  /* 52 */  {9,   4503599627370495U   },/*8191.3.2731.5.53.157.1613*/
  /* 53 */  {71,  9007199254740991U   },/*20394401.6361.69431*/
  /* 54 */  {125, 18014398509481983U  },/*7.262657.73.19.3.3.3.3.87211*/
  /* 55 */  {71,  36028797018963967U  },/*881.23.89.31.3191.201961*/
  /* 56 */  {149, 72057594037927935U  },/*127.43.3.29.113.5.17.15790321*/
  /* 57 */  {45,  144115188075855871U },/*7.524287.1212847.32377*/
  /* 58 */  {99,  288230376151711743U },/*2089.233.1103.3.59.3033169*/
  /* 59 */  {123, 576460752303423487U },/*179951.3203431780337*/
  /* 60 */  { 3,  1152921504606846975U},/*7.151.31.3.3.11.331.13.5.5.41.61.1321*/
  /* 61 */  {39,  2305843009213693951U   },/*2305843009213693951*/
  /* 62 */  {105, 4611686018427387903U   },/*2147483647.3.715827883*/
  /* 63,*/  {3 ,  9223372036854775807U   },/*337.7.7.73.127.649657.92737*/
  /* 64 */  {27,  18446744073709551615U  } /*3.5.17.257.65537.641.6700417*/
#endif
};


/* for checking primitivity */
static coslong gcd ( coslong n, coslong m ) {
  if (m == 0) return(n); 
  else return( gcd( m, n%m ) );
}




static coslong computeReverse ( int n, coslong a ) 
{
  int i;
  coslong a_rev = 1;
  for ( i=n-1; i>0; i-- )
     if ( a & (ONE<<i) ) a_rev |= (ONE << (n-i));     
  return a_rev;
}


static coslong multmod( int n, coslong a, coslong b, coslong p ) 
{
  coslong t, rslt;
  coslong top_bit;
  int i;

  rslt = 0;
  t = a; /* t is a*x^i */
  top_bit = ONE << (n-1);

  for ( i=0; i<n; i++ ) {
     if (b & ONE) rslt ^= t;
     if (t & top_bit)
        t = ( (t & ~top_bit) << 1 ) ^ p;
     else
        t = t << 1;
     b = b >> 1;
  }
  return rslt;
}

static coslong powmod( int n, coslong a, coslong power, coslong p ) 
{
  coslong t = a;
  coslong rslt = ONE;
  while ( power != 0 ) {
    if ( power & ONE ) rslt = multmod( n, t, rslt, p );
    t = multmod( n, t, t, p );
    power = power >> 1;
  }
  return rslt;
}

static coslong minpoly( int n, coslong necklace, coslong p ) 
{
  coslong root, rslt = 0;
  coslong f[ MAX ];
  int i, j;

  f[0] = ONE;
  for (i=1; i<n; i++ ) f[i] = 0;

  root = powmod( n, TWO, necklace, p ); /* '2' is monomial x */
  for (i=1; i<=n; i++ ) {
    if (i != 1)
      root = multmod( n, root, root, p );

    for (j=n-1; j>=1; j--)
      f[j] = f[j-1] ^ multmod( n, f[j], root, p );
    f[0] = multmod( n, f[j], root, p );
  }

  for (i=0; i<n; i++ )
    if (f[i] == ONE)
      rslt |= ONE << i;
    else if (f[i] != 0)
      fprintf( stderr, "Ahh!" );

  return rslt;
}

static coslong toInt (int n, int *b) 
{
  coslong x;
  int i;
  x = 0;
  for (i=1; i<=n; ++i) {
    x = TWO*x;
    if (b[i] == 1) ++x;
  }
  return x;
}


static void OutputPolys(int n, coslong m, coslong gcdResult, coslong *i, coslong imax, coslong *out, int *degrees) 
{    
  if (*i < imax) {
    out[*i] = (ONE<<n) | m; /*computeReverse(n, m); */
    degrees[*i] = n; 
  }
  (*i)++;
} 

static void PrintIt(coslong p, int n, int *b, int primitiveonly, coslong *i, coslong imax, coslong *out, int *degrees ) 
{
  coslong necklace;
  coslong m, m_rev;
  coslong gcdResult;
  static int count = 0;
  
  if (p != n) return; 

  necklace = toInt(n, b);
  m = minpoly( n, necklace, poly_table[n].poly );

  /* check if it's primitive */
  gcdResult = gcd( necklace, poly_table[n].pow2m1 );
  m_rev = computeReverse( n, m );

  if ( (gcdResult == 1) || !primitiveonly ) { /* it's primitive or we want all irreducables */ 
    OutputPolys( n, m, gcdResult,  i, imax, out, degrees );
    if ( m != m_rev ) { 
       OutputPolys( n, m_rev, gcdResult, i, imax, out, degrees );
    } 
  }
}

static void gen(int t, int p, int c, int *b, int n, int primitiveonly, coslong *i, coslong imax, coslong *out, int *degrees) 
{
  if (imax > 0 && (*i) >= imax) return;
  if (t > n) PrintIt( p, n, b, primitiveonly, i, imax, out, degrees );
  else {
    if (b[t-c] == 0) {
      if (b[t-p] == 0) {
        b[t] = 0;  gen( t+1, p, t, b, n, primitiveonly, i, imax, out, degrees );
      }
      b[t] = 1;
      if (b[t-p] == 1) gen( t+1, p, c, b, n, primitiveonly, i, imax, out, degrees );
      else gen( t+1, t, c, b, n, primitiveonly, i, imax, out, degrees );
    } else {
      b[t] = 0;  
      gen( t+1, p, c, b, n, primitiveonly, i, imax, out, degrees );
    }
  }
}


/* return the first k primitive polynomials and their degrees */ 
void get_primitive_polynomials(int k, int *degrees, int *polys)
{
  int i, n; 
  int b[MAX+1];   /* the necklace */

  assert(sizeof(coslong) == sizeof(int)); 
  i=1; 
  n=1; /* current degree */   
  polys[0]=1; 
  degrees[0]=0; 
  b[0] = b[1] = 0; 
  while (i < k) {
    gen( 2, 1, 1, b, n, 1, &i, k, polys, degrees);
    n++; 
  }
    
}

#ifdef POLYTEST
int coscomp(const void *a, const void *b)
{
  return *((int *) a) > *((int *) b); 
}

int main (int argc, char *argv[] ) {
  int imax=100; 
  int out[imax]; 
  int degree[imax]; 
  int i; 

  get_primitive_polynomials(imax, degree, out); 

  qsort(out, imax, sizeof(coslong),  &coscomp);  
  for (i=0; i < imax; i++) {
    printf("%d %7d %7d\n", i, out[i], degree[i]); 
  }
  return(0);
} 
#endif

