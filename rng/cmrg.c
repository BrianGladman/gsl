#include <config.h>
#include <stdlib.h>

/*  From:
    P. L'Ecuyer, "Combined Multiple Recursive Random Number Generators,"
    to appear in Operations Research, 1996.
    (Preprint obtained as file compmrg.ps from L'Ecuyer's web page.)  */

static const int m1 = 2147483647, m2 = 2145483479;

static const int a12 =   63308, q12 = 33921, r12 = 12979,
                 a13 = -183326, q13 = 11714, r13 =  2883,
                 a21 =   86098, q21 = 24919, r21 =  7417,
                 a23 = -539608, q23 =  3976, r23 =  2071 ;

#define gsl_ran_cmrg_RANDMAX 2147483647 /* m1 */

typedef struct {
    long x10,x11,x12;           /* first component */
    long x20,x21,x22;           /* second component */
} gsl_ran_cmrg_randomState;

#define POSITIVE(x,m) if (x<0) x += m

unsigned long cmrg_get (void * vstate)
{
    int h,p12,p13,p21,p23;
    gsl_ran_cmrg_randomState * state = (gsl_ran_cmrg_randomState *) vstate;

    /* Component 1 */
    h = state->x10 / q13 ; p13 = -a13 * (state->x10 - h * q13) - h * r13;
    h = state->x11 / q12 ; p12 =  a12 * (state->x11 - h * q12) - h * r12;
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

#define LCG(n) ((n)*8121+28411)%134456

void cmrg_set(void * vstate, unsigned int s)
{
    /* An entirely adhoc way of seeding! This does **not** come
       from L'Ecuyer et al */
    gsl_ran_cmrg_randomState * state = (gsl_ran_cmrg_randomState *) vstate;

    if (s==0) s=1;

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

static gsl_ran_cmrg_randomState state = {
    511515612L, 1645048169L, 1860274777L,
     55882945L, 1225790668L, 2055528708L
};
#include "cmrg-state.c"

    
    
    
