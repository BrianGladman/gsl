/* Author: G. Jungman
 * RCS:    $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_chebyshev.h"
#include "gsl_sf_psi.h"


/*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*/

/* Chebyshev fits from SLATEC code for psi(x)

 Series for PSI        on the interval  0.	   to  1.00000D+00
				       with weighted error   2.03E-17
					log weighted error  16.69
			      significant figures required  16.39
				   decimal places required  17.37

 Series for APSI       on the interval  0.	   to  2.50000D-01
				       with weighted error   5.54E-17
					log weighted error  16.26
			      significant figures required  14.42
				   decimal places required  16.86

*/

static double psics_data[23] = {
  -.038057080835217922,
   .491415393029387130, 
  -.056815747821244730,
   .008357821225914313,
  -.001333232857994342,
   .000220313287069308,
  -.000037040238178456,
   .000006283793654854,
  -.000001071263908506,
   .000000183128394654,
  -.000000031353509361,
   .000000005372808776,
  -.000000000921168141,
   .000000000157981265,
  -.000000000027098646,
   .000000000004648722,
  -.000000000000797527,
   .000000000000136827,
  -.000000000000023475,
   .000000000000004027,
  -.000000000000000691,
   .000000000000000118,
  -.000000000000000020
};    
static double apsics_data[16] = {    
  -.0204749044678185,
  -.0101801271534859,
   .0000559718725387,
  -.0000012917176570,
   .0000000572858606,
  -.0000000038213539,
   .0000000003397434,
  -.0000000000374838,
   .0000000000048990,
  -.0000000000007344,
   .0000000000001233,
  -.0000000000000228,
   .0000000000000045,
  -.0000000000000009,
   .0000000000000002,
  -.0000000000000000 
};    
static struct gsl_sf_cheb_series psi_cs = {
  psics_data,
  22,
  -1, 1,
  (double *)0,
  (double *)0
};
static struct gsl_sf_cheb_series apsi_cs = {
  apsics_data,
  15,
  -1, 1,
  (double *)0,
  (double *)0
};


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

/* simple implementation for integer argument */
int gsl_sf_psi_int_impl(const int n, double * result)
{
  if(n == 1) {
    *result = -M_EULER;
    return GSL_SUCCESS;
  }
  else if(n <= 0) {
    *result = 0.;
    return GSL_EDOM;
  }
  else {
    int k;
    double ans;
    ans = -M_EULER;
    for(k=1; k<n; k++) {
      ans += 1.0/k;
    }
    *result = ans;
    return GSL_SUCCESS;
  }
}



int gsl_sf_psi_impl(const double x, double * result)
{
  double y = fabs(x);
  double xbig  = 1.0  / GSL_SQRT_MACH_EPS;    /* XBIG  = 1.0/SQRT(R1MACH(3)) */
  double dxrel = 10.0 * GSL_SQRT_MACH_EPS;    /* DXREL = SQRT (R1MACH(4))    */
  
  if(y >= 2.0) {
    double aux = 0.;
    if(y < xbig) aux = gsl_sf_cheb_eval(8./(y*y)-1., &apsi_cs);
    if(x < 0.0) {
      /* *result = log(y) - 0.5/x + aux - M_PI * cot(M_PI*x); */
      *result = log(y) - 0.5/x + aux - M_PI * cos(M_PI*x)/sin(M_PI*x);
    }
    else {
      *result = log(y) - 0.5/x + aux;
    }
    return GSL_SUCCESS;
  }
  else { /* y < 2.0 */
    if(x == 0.0) {
      *result = 0.0; /* FIXME: should be Inf */
      return GSL_EDOM;
    }
    else {
      double ans;
      int n = x;
      if(x < 0.0) --n;
      y = x - n;
      --n;
      ans = gsl_sf_cheb_eval(2.*y-1., &psi_cs);
      if(n == 0) {
	*result = ans;
	return GSL_SUCCESS;
      }

      n = -n;

      if(x < 0.0 && x+n-2 == 0.) {
      	/* x is a negative integer */
	*result = 0.; /* FIXME: should be Inf */
	return GSL_EDOM;
      }
      else {
	int i;
	for(i=0; i<n; i++) {
          ans -= 1./(x + i);
      	}
	*result = ans;
	if(x < -0.5 && fabs((x-(int)(x-0.5))/x) < dxrel) {
      	  /* loss of precision: x near a negative integer */
	  return GSL_ELOSS;
      	}
	else {
	  return GSL_SUCCESS;
	}
      }
    }
  }
}


/*-*-*-*-*-*-*-*s = 0;
  s += ( frac_diff(gsl_sf_psi(5000.0), 8.517093188082904107) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_psi(5000.0)");
  status += s;-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_psi_int_e(const int n, double * result)
{
  int status = gsl_sf_psi_int_impl(n, result);
  if(status != GSL_SUCCESS){
    GSL_ERROR("gsl_sf_psi_int_e", status);
  }
  return status;
}

int gsl_sf_psi_e(const double x, double * result)
{
  int status = gsl_sf_psi_impl(x, result);
  if(status != GSL_SUCCESS){
    GSL_ERROR("gsl_sf_psi_e", status);
  }
  return status;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*-*/

double gsl_sf_psi_int(const int n)
{
  double y;
  int status = gsl_sf_psi_int_impl(n, &y);
  if(status != GSL_SUCCESS){
    GSL_WARNING("gsl_sf_psi_int", status);
  }
  return y;
}

double gsl_sf_psi(const double x)
{
  double y;
  int status = gsl_sf_psi_impl(x, &y);
  if(status != GSL_SUCCESS){
    GSL_WARNING("gsl_sf_psi", status);
  }
  return y;
}
