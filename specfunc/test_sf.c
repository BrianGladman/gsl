/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <stdio.h>
#include <stdlib.h>
#include <gsl_sf.h>


double frac_diff(double x1, double x2)
{
  return fabs((x1-x2)/(x1+x2));
}


int check_airy(void)
{
  int status = 0;
  int s;

  s = 0;
  s += ( frac_diff(gsl_sf_airy_Ai(0.), 1.) > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Ai(0.), 1.) > 1.e-14 );
  gsl_test(s, "gsl_sf_airy_Ai");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_airy_Bi(0.), 1.) > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Bi(0.), 1.) > 1.e-14 );
  gsl_test(s, "gsl_sf_airy_Bi");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_airy_Bi_scaled(0.), 1.) > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Bi_scaled(0.), 1.) > 1.e-14 );
  gsl_test(s, "gsl_sf_airy_Bi_scaled");
  status += s;

  return status;
}


int check_bessel(void)
{
  int status = 0;
  int s;

  return status;
}


int main(int argc, char * argv[])
{

  gsl_test(check_airy(),       "Airy functions")
  gsl_test(check_bessel(),     "Bessel functions")  

  exit(0);  
}
