#include <math.h>
#include "gsl_specfunc.h"


void complex_log(double zr, double zi, double * lnr, double * theta)
{
  if(zr != 0.0 || zi != 0.0) {
    double r2 = zr*zr + zi*zi;
    *lnr = 0.5*log(r2);
    *theta = atan2(zi, zr);
  }
  else {
    char buff[50];
    sprintf(buff,"complex_log: z=0.0");
    push_error(buff, Error_Domain_);
  }
}
