/* $Id$ */

#include <stdio.h>
#include <stdlib.h>
#include "gsl_special.h"

const double LOGBASETENOFTWO=0.301029995663981;
const double SQRTONEHALF=0.707106781186547;
const double ONEOVERLNTEN=0.434294481903252;
double sig2lamp(double sig)
{
    return LOGBASETENOFTWO - ONEOVERLNTEN*GSL_log_erfc(sig*SQRTONEHALF);
}

int
main(int argc, char **argv)
{
    double x;
    double xo,xf,dx;

    xo=xf=dx=1.0;

    if (argc>1) xo=atof(argv[1]);
    if (argc>2) xf=atof(argv[2]);
    if (argc>3) dx=atof(argv[3]);

    for (x=xo; x<=xf; x+=dx) {
        printf("%g %g\n",x,sig2lamp(x));
    }
}
