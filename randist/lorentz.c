


/* $Id$ */
#include <math.h>
#include "gsl_ran.h"
#include "gsl_randist.h"

double 
gsl_ran_lorentzian (double mu)
{
  double u;
  do
    {
      u = gsl_ran_uniform ();
    }
  while (u == 0.5);

  return mu * tan (M_PI * u);
}

/*  The Lorentzian probability distribution is 

   1
   p(x) =  -----------------
   pi (1 + (x/mu)^2) 

   It is also known as the Cauchy distribution */
