#include <config.h>
#include <stdlib.h>
#include <gsl_rng.h>

/*  This is a combined multiple recursive generator. The sequence is,

         z_n = (x_n - y_n) mod m1

    where,
         
         x_n = (a_{11} x_{n-1} + a_{12} x_{n-2} + a_{13{ x_{n-3}) mod m1
         y_n = (a_{21} y_{n-1} + a_{22{ y_{n-2} + a_{23} y_{n-3}) mod m2
    
         a_{11} = 0, a_{12} = 63308, a_{13} = -183326
         a_{21} = 86098, a_{22} = 0, a_{23} = -539608
         
         m1 = 2^31 - 1 = 2147483647
         m2 = 2^31 - 2000169 = 2145483479

    We seed the generator with

        x0 = (A * seed + B) mod M
        x1 = (A * x0 + B) mod M
        x2 = (A * x1 + B) mod M
        y0 = (A * x2 + B) mod M
        y1 = (A * y0 + B) mod M
        y2 = (A * y1 + B) mod M

    and then use 7 iterations of the generator to "warm up" the
    internal state.

    For checking the theoretical value of z_{10008} is 1477798470.
    The subscript 10008 means (1) seed the generator with s=1, (2) do
    seven warm-up iterations, (3) then do 10000 actual iterations.

    The period of this generator is about 2^{205}.

    From: P. L'Ecuyer, "Combined Multiple Recursive Random Number
    Generators," Operations Research, 44, 5 (1996), 816--822.

    This is available on the net from L'Ecuyer's home page,

    http://www.iro.umontreal.ca/~lecuyer/myftp/papers/combmrg.ps
    ftp://ftp.iro.umontreal.ca/pub/simulation/lecuyer/papers/combmrg.ps */

unsigned long int cmrg_get (void * vstate);
void cmrg_set (void * state, unsigned int s);

static const int m1 = 2147483647, m2 = 2145483479;

static const int a12 =   63308, q12 = 33921, r12 = 12979,
                 a13 = -183326, q13 = 11714, r13 =  2883,
                 a21 =   86098, q21 = 24919, r21 =  7417,
                 a23 = -539608, q23 =  3976, r23 =  2071 ;

typedef struct {
  long int x10, x11, x12;           /* first component */
  long int x20, x21, x22;           /* second component */
} cmrg_state_t;

unsigned long int cmrg_get (void * vstate)
{
  int h,p12,p13,p21,p23;
  cmrg_state_t * state = (cmrg_state_t *) vstate;
  
  /* Component 1 */
  h = state->x10 / q13 ; p13 = -a13 * (state->x10 - h * q13) - h * r13;
  h = state->x11 / q12 ; p12 =  a12 * (state->x11 - h * q12) - h * r12;
#define POSITIVE(x,m) if (x<0) x += m
  POSITIVE(p13,m1);
  POSITIVE(p12,m1);
  state->x10 = state->x11;
  state->x11 = state->x12;
  state->x12 = p12 - p13;
  POSITIVE(state->x12,m1);
  
  /* Component 2 */
  h = state->x20 / q23 ; p23 = -a23 * (state->x20 - h * q23)  - h * r23;
  h = state->x22 / q21 ; p21 =  a21 * (state->x22 - h * q21)  - h * r21;
  POSITIVE(p23,m2);
  POSITIVE(p21,m2);
  state->x20 = state->x21;
  state->x21 = state->x22;
  state->x22 = p21 - p23;
  POSITIVE(state->x22,m2);
  
  /* Combination */
  if(state->x12 < state->x22)
    return (state->x12 - state->x22 + m1);
  else
    return (state->x12 - state->x22);
}

void cmrg_set(void * vstate, unsigned int s)
{
  /* An entirely adhoc way of seeding! This does **not** come from
     L'Ecuyer et al */
  
  cmrg_state_t * state = (cmrg_state_t *) vstate;
  
  if (s == 0) s = 1;
  
#define LCG(n) ((n)*8121+28411)%134456
  state->x10 = LCG(s);
  state->x11 = LCG(state->x10);
  state->x12 = LCG(state->x11);
  state->x20 = LCG(state->x12);
  state->x21 = LCG(state->x20);
  state->x22 = LCG(state->x21);
  
  /* "warm it up" */
  cmrg_get (state);
  cmrg_get (state);
  cmrg_get (state);
  cmrg_get (state);
  cmrg_get (state);
  cmrg_get (state);
  cmrg_get (state);
}

static const gsl_rng_type cmrg_type = { "cmrg",  /* name */
					2147483647,  /* RAND_MAX */
					sizeof(cmrg_state_t), 
					&cmrg_set, 
					&cmrg_get } ;

const gsl_rng_type * gsl_rng_cmrg = &cmrg_type ;



    
    
    
