/* $Id$ */

#include <stdio.h>
#include "gsl_sf.h"

int
main(int argc, char **argv)
{
    double x;

    for (x=0; x<30; x+=0.1) {
        printf("%g %g\n",x,GSL_Q(x));
    }
    for (x=0; x<30; x+=0.1) {
        printf("%g %g\n",x,0.5*GSL_erfc(x*0.707106781186547));
    }
    
    /*
    for (x=2.0; x<200.0; x+=0.1) {
        printf("%g %g\n",x,GSL_log_erfc_asymptotic(x)-GSL_log_erfc(x));
    }
    for (x=-3.0; x<=3.0; x+=0.01) {
        printf("%g %g\n",x,GSL_erfc(x)+GSL_erf(x)-1.0);
    }
    printf("#k\n");
    for (x=0; x<10; x+=0.1) {
        printf("%g %g\n",x,1.0-GSL_erf(x));
    }
    */
}
