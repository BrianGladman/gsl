/* Author: G. Jungman
 * RCS:    $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_chebyshev.h"
#include "gsl_sf_psi.h"


/*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*/


/* Chebyshev fit for f(y) = Re(Psi(1+Iy)) + M_EULER - y^2/(1+y^2) - y^2/(2(4+y^2))
 * 1 < y < 10
 *   ==>
 * y(x) = (9x + 11)/2,  -1 < x < 1
 * x(y) = (2y - 11)/9
 *
 * g(x) := f(y(x))
 */
static double r1py_data[] = {
   1.59888328244976954803168395603,
   0.67905625353213463845115658455,
  -0.068485802980122530009506482524,
  -0.0057881841830958667920088311823,
   0.0085112581671086159804198556475,
  -0.0040426561346996934343345564091,
   0.00135232840615940260177846295622,
  -0.000311646563930660566674525382102,
   0.0000185075637852491354372191392113,
   0.0000283487054275298502964921455903,
  -0.0000194875360145745355675419596539,
   8.0709788710834469408621587335e-6,
  -2.29835643213405180370603465611e-6,
   3.05066295996047498438559626587e-7,
   1.30422386324183646107742848462e-7,
  -1.23086571810489505894646902083e-7,
   5.7710855710682427240667414345e-8,
  -1.82755593424509639660926363536e-8,
   3.10204713006265894207595189301e-9,
   6.8989327480593812470039430640e-10,
  -8.7182290258923059852334818997e-10,
   4.4069147710243611798213548777e-10,
  -1.47273110991985359634672002769e-10,
   2.75896825232626447488258442482e-11,
   4.1871826756975856411554363568e-12,
  -6.5673460487260087541400767340e-12,
   3.4487900886723214020103638000e-12,
  -1.18072514174486906079737940779e-12,
   2.37983143439695892587093155740e-13,
   2.16636304108188318242594658208e-15
};
static struct gsl_sf_cheb_series r1py_cs = {
  r1py_data,
  29,
  -1,1,
  (double *)0,
  (double *)0
};


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
    if(y < xbig) aux = gsl_sf_cheb_eval(&apsi_cs, 8.0/(y*y)-1.0);
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
      ans = gsl_sf_cheb_eval(&psi_cs, 2.0*y-1.0);
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


int
gsl_sf_psi_1piy_impl(const double y, double * result)
{
  double ay = fabs(y);

  if(ay > 1000.0) {
    /* [Abramowitz+Stegun, 6.3.19] */
    double yi2 = 1.0/(ay*ay);
    double ln  = log(y);
    double sum = yi2 * (1.0/12.0 + 1.0/120.0 * yi2 + 1.0/252.0 * yi2*yi2);
    *result = ln + sum;
    return GSL_SUCCESS;
  }
  else if(ay > 10.0) {
    /* [Abramowitz+Stegun, 6.3.19] */
    double yi2 = 1.0/(ay*ay);
    double ln  = log(y);
    double sum = yi2 * (1.0/12.0 +
                   yi2 * (1.0/120.0 +
		     yi2 * (1.0/252.0 +
                       yi2 * (1.0/240.0 +
		         yi2 * (1.0/132.0 + 691.0/32760.0 * yi2)))));
    *result = ln + sum;
    return GSL_SUCCESS;
  }
  else if(ay > 1.0){
    double y2 = ay*ay;
    double x = (2.0*ay - 11.0)/9.0;
    double c = gsl_sf_cheb_eval(&r1py_cs, x);
    *result = c - M_EULER + y2*(1.0/(1.0+y2) + 0.5/(4.0+y2));
    return GSL_SUCCESS;
  }
  else {
    /* [Abramowitz+Stegun, 6.3.17]
     *
     * -M_EULER + y^2 Sum[1/n 1/(n^2 + y^2), {n,1,M}]
     *   +     Sum[1/n^3, {n,M+1,Infinity}]
     *   - y^2 Sum[1/n^5, {n,M+1,Infinity}]
     *   + y^4 Sum[1/n^7, {n,M+1,Infinity}]
     *   - y^6 Sum[1/n^9, {n,M+1,Infinity}]
     *   + O(y^8)
     *
     * We take M=50 for at least 15 digit precision.
     */
    const int M = 50;
    const double y2 = y*y;
    const double c0 = 0.00019603999466879846570;
    const double c2 = 3.8426659205114376860e-08;
    const double c4 = 1.0041592839497643554e-11;
    const double c6 = 2.9516743763500191289e-15;
    const double p  = c0 + y2 *(-c2 + y2*(c4 - y2*c6));
    double sum = 0.0;
    
    int n;
    for(n=1; n<=M; n++) {
      sum += 1.0/(n * (n*n + y*y));
    }
    
    *result = -M_EULER + y2 * (sum + p);
    return GSL_SUCCESS;
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

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

int gsl_sf_psi_1piy_e(const double x, double * result)
{
  int status = gsl_sf_psi_1piy_impl(x, result);
  if(status != GSL_SUCCESS){
    GSL_ERROR("gsl_sf_psi_1piy_e", status);
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

double gsl_sf_psi_1piy(const double x)
{
  double y;
  int status = gsl_sf_psi_1piy_impl(x, &y);
  if(status != GSL_SUCCESS){
    GSL_WARNING("gsl_sf_psi_1piy", status);
  }
  return y;
}
