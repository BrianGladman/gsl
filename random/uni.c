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

#include <stdlib.h>
#include "uni.h"

static const int MDIG=16;                /* Machine digits in int */
static const int m1 = 32767;             /* 2^(MDIG-1) - 1 */
static const int m2 = 256;               /* 2^(MDIG/2) */

#define gsl_ran_uni_RANDMAX 32767 /*m1*/

typedef struct {
    int i,j;
    unsigned long m[17];
} gsl_ran_uni_randomState;

static void
gsl_ran_uni_printState_p(gsl_ran_uni_randomState *s)
{
    int n;
    printf("%d, %d,  {\n",s->i,s->j);
    for (n=0; n<16; ++n) {
	printf("%lu,%c",s->m[n],(n%5==4 ? '\n' : ' '));
    }
    printf("%lu }\n",s->m[16]);
}


inline unsigned long gsl_ran_uni_random_wstate(void *vState)
{
    long k;                     /* important k not be unsigned */
    gsl_ran_uni_randomState *theState;
    theState = (gsl_ran_uni_randomState *)vState;

    k = theState->m[theState->i] - theState->m[theState->j];
    if (k < 0) k += m1;
    theState->m[theState->j] = k;
    
    if (--(theState->i) == -1) theState->i=16;
    if (--(theState->j) == -1) theState->j=16;

    return k;

}

void gsl_ran_uni_seed_wstate(void *vState, int jd)
{
    int i,jseed,k0,k1,j0,j1;
    gsl_ran_uni_randomState *theState;
    theState = (gsl_ran_uni_randomState *)vState;
    
    /* For this routine, the seeding is very elaborate! */
    /* A flaw in this approach is that seeds 1,2 give exactly the
       same random number sequence!  */

    jd = (jd > 0 ? jd : -jd);       /* absolute value */
    jd = 2*jd+1;                    /* enforce seed be odd */
    jseed = (jd < m1 ? jd : m1);    /* seed should be less than m1 */
                                    
    k0 = 9069%m2;
    k1 = 9069/m2;
    j0 = jseed%m2;
    j1 = jseed/m2;

    for (i=0; i<17; ++i) {
        jseed = j0*k0;
        j1 = (jseed/m2 + j0*k1 + j1*k0) % (m2/2);
        j0 = jseed%m2;
        theState->m[i] = j0+m2*j1;
    }
    theState->i=4;
    theState->j=16;

    return;
}


static gsl_ran_uni_randomState state = {
    4, 16,  {
	27207, 30011, 31519, 10547, 951,
	6635, 10767, 30051, 1063, 6555,
	6143, 5267, 23447, 9291, 13551,
	14019, 31239 }
/* The numbers below were provided in the version that came
 * over the net.  I have used the numbers above following the
 * convention that the default initializer is the same as you
 * would get if you used seed=1
    4, 16, {
	30788, 23052,  2053, 19346, 10646, 19427, 23975,
	19049, 10949, 19693, 29746, 26748, 2796,  23890,
	29168, 31924, 16499 }
 */
};
#include "uni-state.c"


    
    
    
