/**
  This is a lagged Fibonacci generator which supposedly excellent
  statistical properties (I do not concur)

  I got it from the net and translated into C.

* ======================================================================
* NIST Guide to Available Math Software.
* Fullsource for module UNI from package CMLIB.
* Retrieved from CAMSUN on Tue Oct  8 14:04:10 1996.
* ======================================================================

C***BEGIN PROLOGUE  UNI
C***DATE WRITTEN   810915
C***REVISION DATE  830805
C***CATEGORY NO.  L6A21
C***KEYWORDS  RANDOM NUMBERS, UNIFORM RANDOM NUMBERS
C***AUTHOR    BLUE, JAMES, SCIENTIFIC COMPUTING DIVISION, NBS
C             KAHANER, DAVID, SCIENTIFIC COMPUTING DIVISION, NBS
C             MARSAGLIA, GEORGE, COMPUTER SCIENCE DEPT., WASH STATE UNIV
C
C***PURPOSE  THIS ROUTINE GENERATES QUASI UNIFORM RANDOM NUMBERS ON [0,1
C             AND CAN BE USED ON ANY COMPUTER WITH WHICH ALLOWS INTEGERS
C             AT LEAST AS LARGE AS 32767.
C***DESCRIPTION
C
C       THIS ROUTINE GENERATES QUASI UNIFORM RANDOM NUMBERS ON THE INTER
C       [0,1).  IT CAN BE USED WITH ANY COMPUTER WHICH ALLOWS
C       INTEGERS AT LEAST AS LARGE AS 32767.
C
C
C   USE
C       FIRST TIME....
C                   Z = UNI(JD)
C                     HERE JD IS ANY  N O N - Z E R O  INTEGER.
C                     THIS CAUSES INITIALIZATION OF THE PROGRAM
C                     AND THE FIRST RANDOM NUMBER TO BE RETURNED AS Z.
C       SUBSEQUENT TIMES...
C                   Z = UNI(0)
C                     CAUSES THE NEXT RANDOM NUMBER TO BE RETURNED AS Z.
C
C
C..................................................................
C   NOTE: USERS WHO WISH TO TRANSPORT THIS PROGRAM FROM ONE COMPUTER
C         TO ANOTHER SHOULD READ THE FOLLOWING INFORMATION.....
C
C   MACHINE DEPENDENCIES...
C      MDIG = A LOWER BOUND ON THE NUMBER OF BINARY DIGITS AVAILABLE
C              FOR REPRESENTING INTEGERS, INCLUDING THE SIGN BIT.
C              THIS VALUE MUST BE AT LEAST 16, BUT MAY BE INCREASED
C              IN LINE WITH REMARK A BELOW.
C
C   REMARKS...
C     A. THIS PROGRAM CAN BE USED IN TWO WAYS:
C        (1) TO OBTAIN REPEATABLE RESULTS ON DIFFERENT COMPUTERS,
C            SET 'MDIG' TO THE SMALLEST OF ITS VALUES ON EACH, OR,
C        (2) TO ALLOW THE LONGEST SEQUENCE OF RANDOM NUMBERS TO BE
C            GENERATED WITHOUT CYCLING (REPEATING) SET 'MDIG' TO THE
C            LARGEST POSSIBLE VALUE.
C     B. THE SEQUENCE OF NUMBERS GENERATED DEPENDS ON THE INITIAL
C          INPUT 'JD' AS WELL AS THE VALUE OF 'MDIG'.
C          IF MDIG=16 ONE SHOULD FIND THAT
   Editors Note: set the seed using 152 in order to get uni(305)
   -jt
C            THE FIRST EVALUATION
C              Z=UNI(305) GIVES Z=.027832881...
C            THE SECOND EVALUATION
C              Z=UNI(0) GIVES   Z=.56102176...
C            THE THIRD EVALUATION
C              Z=UNI(0) GIVES   Z=.41456343...
C            THE THOUSANDTH EVALUATION
C              Z=UNI(0) GIVES   Z=.19797357...
C
C***REFERENCES  MARSAGLIA G., "COMMENTS ON THE PERFECT UNIFORM RANDOM
C                 NUMBER GENERATOR", UNPUBLISHED NOTES, WASH S. U.
C***ROUTINES CALLED  I1MACH,XERROR
C***END PROLOGUE  UNI

  **/

#include <config.h>
#include <stdlib.h>
#include <gsl_rng.h>

unsigned long int uni_get (void * vstate);
void uni_set (void * state, unsigned int s);
void uni_set_with_state (void * vstate, const void * vinit_state,
			  unsigned int s);

static const int MDIG=16;                /* Machine digits in int */
static const int m1 = 32767;             /* 2^(MDIG-1) - 1 */
static const int m2 = 256;               /* 2^(MDIG/2) */

typedef struct {
    int i,j;
    unsigned long m[17];
} uni_state_t;

inline unsigned long uni_get (void * vstate)
{
    uni_state_t * state = (uni_state_t *) vstate;
    const int i = state->i ;
    const int j = state->j ;

    /* important k not be unsigned */
    long k = state->m[i] - state->m[j];

    if (k < 0) k += m1;
    state->m[j] = k;
    
    if (i == 0) 
      {
	state->i = 16 ;
      } 
    else
      {
	(state->i)-- ;
      }

    if (j == 0) 
      {
	state->j = 16 ;
      } 
    else
      {
	(state->j)-- ;
      }

    return k;
}

static const uni_state_t init_state = {
  4, 16,  { 27207, 30011, 31519, 10547, 951,
	    6635, 10767, 30051, 1063, 6555,
	    6143, 5267, 23447, 9291, 13551,
	    14019, 31239 }
  /* The numbers below were provided in the version that came
     over the net.  I have used the numbers above following the
     convention that the default initializer is the same as you
     would get if you used seed=1
   4, 16, { 30788, 23052,  2053, 19346, 10646, 19427, 23975,
            19049, 10949, 19693, 29746, 26748, 2796,  23890,
	    29168, 31924, 16499 } */
};

void uni_set(void * state, unsigned int s)
{
  uni_set_with_state (state, &init_state, s) ;
}

void uni_set_with_state (void * vstate, const void * vinit_state, 
			 unsigned int s)
{
  int i,seed,k0,k1,j0,j1;
  
  uni_state_t * state = (uni_state_t *) vstate;
  
  *state = *(const uni_state_t *) vinit_state ;
    
  /* For this routine, the seeding is very elaborate! */
  /* A flaw in this approach is that seeds 1,2 give exactly the
     same random number sequence!  */
  
  s = 2*s+1;                    /* enforce seed be odd */
  seed = (s < m1 ? s : m1);    /* seed should be less than m1 */
  
  k0 = 9069%m2;
  k1 = 9069/m2;
  j0 = seed%m2;
  j1 = seed/m2;
  
  for (i=0; i<17; ++i) {
    seed = j0*k0;
    j1 = (seed/m2 + j0*k1 + j1*k0) % (m2/2);
    j0 = seed%m2;
    state->m[i] = j0+m2*j1;
  }
  state->i=4;
  state->j=16;
  
  return;
}

static const gsl_rng_type uni_type = { "gsl-uni",  /* name */
					32767,  /* RAND_MAX */
					sizeof(uni_state_t), 
					&uni_set, 
					&uni_get } ;

const gsl_rng_type * gsl_rng_uni (void) { return &uni_type ; }






    
    
    
