/* Gamma distribution of order a>0 is defined by:
 *
 *   F(x) = \frac{1}{\Gamma(a)}\int_0^x t^{a-1} e^{-t} dt
 *
 *   for x>0.  If X and Y are independent gamma-distributed
 *   random variables of order a and b, then X+Y has gamma
 *   distribution of order a+b.
 *
 *   The algorithms below are from Knuth, vol 2, 2nd ed, p. 129.
 */

#include <math.h>
#include "gsl_ran.h"
#include "gsl_randist.h"

#define GSL_PI 3.14159265358979323844
#define GSL_E  2.71828182845904523536
#define GSL_LOGINFINITY 300.0

static double gsl_ran_gammalarge(double a)
{
    /* Works only if a>1, and is most efficient if a is large */
    /* This algorithm, reported in Knuth, is attributed to Ahrens.  A
       faster one, we are told, can be found in: J. H. Ahrens and
       U. Dieter, Computing 12 (1974) 223-246.  */
    double sqa,x,y,v;
    sqa = sqrt(2*a-1);
    do {
        do {
            y = tan(GSL_PI*gsl_ran_uniform());
            x = sqa*y + a-1;
        } while (x <= 0);
        v = gsl_ran_uniform();
    } while (v > (1+y*y)*exp((a-1)*log(x/(a-1))-sqa*y));

    return x;
}

static double gsl_ran_gammafrac(double a)
{
    /* This is exercise 16 from Knuth; see page 135,
     * and the solution is on page 551.
     */
    double p,q,x,u,v;
    p = GSL_E/(a+GSL_E);
    do {
        u = gsl_ran_uniform();
        do {
            v = gsl_ran_uniform();
        } while (v==0);
        if (u < p) {
            x = exp((1/a)*log(v));
            q = exp(-x);
        } else {
            x = 1 - log(v);
            q = exp((a-1)*log(x));
        }
    } while( gsl_ran_uniform() >= q );
    return x;
}

double gsl_ran_gamma(double a)
{
    /* assume a>0 */
    int na;
    na = floor(a);
    if (a == na) {
        return gsl_ran_gammaint(na);
    } else if (na == 0) {
        return gsl_ran_gammafrac(a);
    } else {
        return gsl_ran_gammaint(na) + gsl_ran_gammafrac(a-na);
    }
}

double gsl_ran_gammaint(int a)
{
    if (a < 12) {
        int i;
        double prod=1.0;
        for (i=0; i<a; ++i)
            prod *= gsl_ran_uniform();
        if (prod == 0) {
            return GSL_LOGINFINITY;
        } else {
            return -log(prod);
        }
    } else {
        return gsl_ran_gammalarge((double)a);
    }
}
        
        
