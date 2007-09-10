/* integration/integration.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Created: [GJ] Tue Apr 23 21:26:53 EDT 1996
 */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_errno.h>
/* #include "interpolation.h" */
#include <gsl/gsl_integration.h>
#if 0


static char * info_string = 0;

void set_integ_info(const char * mess)
{
  if(info_string != 0) free(info_string);
  info_string = 0;
  if(mess != 0) {
    info_string = (char *)malloc((strlen(mess)+1) * sizeof(char));
    strcpy(info_string, mess);
  }
}


/* Integration function for use with open_romberg(). */
#define FUNC(x) ((*func)(x))
static double midpnt(double (*func)(double), double a, double b, int n)
{
  double x,tnm,sum,del,ddel;
  static double s;
  static int it;
  int j;

  if (n == 1) {
    it=1;
    return (s=(b-a)*FUNC(0.5*(a+b)));
  }
  else {
    tnm=it;
    del=(b-a)/(3.0*tnm);
    ddel=del+del;
    x=a+0.5*del;
    sum=0.0;
    for (j=1;j<=it;j++) {
      sum += FUNC(x);
      x += ddel;
      sum += FUNC(x);
      x += del;
    }
    it *= 3;
    s=(s+(b-a)*sum/tnm)/3.0;
    return s;
  }
}
#undef FUNC


#define JMAX 10
#define JMAXP JMAX+1
#define K 5
double open_romberg(double(*func)(double), double a, double b, double eps)
{
  int j;
  double ss,dss,h[JMAXP+1],s[JMAXP+1];

  h[1]=1.0;
  for (j=1;j<=JMAX;j++) {
    s[j]=midpnt(func,a,b,j);
    if (j >= K) {
      /* local_polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss); */
      interp_poly(&h[j-K]+1,&s[j-K]+1,K,0.0,&ss,&dss);
      if (fabs(dss) < eps * fabs(ss)) return ss;
    }
    s[j+1]=s[j];
    h[j+1]=h[j]/9.0;
  }

  push_error("open_romberg: too many steps", Error_ConvFail_);
  push_generic_error("open_romberg:", info_string);

  return 0.;
}
#undef JMAX
#undef JMAXP
#undef K


double gauss_legendre_10(double (*func)(double), double a, double b)
{
  int j;
  static double x[] = {0., 
                         0.1488743389,
                         0.4333953941,
                         0.6794095682,
                         0.8650633666,
                         0.9739065285};
  static double w[] = {0.,
                         0.2955242247,
                         0.2692667193,
                         0.2190863625,
                         0.1494513491,
                         0.0666713443};
  
  double xm = 0.5 * (b + a);
  double xr = 0.5 * (b - a);
  
  double s = 0.;
  double dx;
  double result;
  
  for(j=1; j<=5; j++){
    dx = xr * x[j];
    s += w[j] * ((*func)(xm+dx) + (*func)(xm-dx));
  }
  result = s * xr;
  
  return result;
}


/************************************************************************
 *                                                                      *
 * Trapezoid rule.                                                      *
 *                                                                      *
 * The original trapzd() function from the Numerical Recipes worked     *
 * by side-effects; it tried to remember the value of the integral      *
 * at the current level of refinement. This is REAL BAD, because if     *
 * you try to use it to do a double integral you will get a surprise.   *
 * You cannot tell which "integral" it is remembering. This stems from  *
 * the stupid fact that there was only one, essentially global,         *
 * variable that was doing this memory job.                             *
 *                                                                      *
 * The solution is simple: pass the current refinement to the function  *
 * so that it can be remembered by an external environment, making it   *
 * easy to avoid confusion.                                             *
 *                                                                      *
 * So the new-method code-fragment for doing an integral looks like:    *
 *                                                                      *
 *       double answer;                                                 *
 *       for(j=1; j<=M+1; j++)                                          *
 *         trapezoid_rule(func, a, b, j, &answer);                      *
 *                                                                      *
 ************************************************************************/
void trapezoid_rule(double(*f)(double), double a, double b, int n, double *s)
{
  double x, tnm, sum, del;
  int it, j;

  if(n==1){
    *s = 0.5 * (b-a) * (f(b) + f(a));
  }
  else {
    for(it=1, j=1; j < n-1; j++)  it <<= 1;
    tnm = (double) it;
    del = (b-a) / tnm;
    x = a + 0.5 * del;
    
    for(sum=0., j=1; j<=it; j++, x+=del) { sum += f(x); }

    *s = 0.5 * (*s + del * sum);
  }
}


/************************************************************************
 *                                                                      *
 * Trapezoid rule.                                                      *
 * This version produces a tracing output.                              *
 *                                                                      *
 ************************************************************************/
#define FUNC(x) ((*func)(x))

void test_trapezoid_rule(double(*func)(double), double a, double b, int n,
                         double *s)
{
  double x, tnm, sum, del;
  int it, j;

  if(n==1){
    printf("t:  a= %g  b= %g   f(a)= %g  f(b)= %g\n",
           a, b, FUNC(a), FUNC(b));
  }

  if(n==1){
    *s = 0.5 * (b-a) * (FUNC(b) + FUNC(a));
    printf("s= %g\n", *s);
  }
  else {
    for(it=1, j=1; j < n-1; j++)  it <<= 1;
    tnm = (double) it;
    del = (b-a) / tnm;
    x = a + 0.5 * del;
    
    for(sum=0., j=1; j<=it; j++, x+=del) sum += FUNC(x);

    *s = 0.5 * (*s + del * sum);

    printf("sum= %g   tnm= %g  del= %g  s= %g\n", sum, tnm, del, *s);
  }
}
#undef FUNC



/* This is fixed to use the non-side-effecting version of the
 * trapezoidal rule, as implemented in trapezoid_rule() above.
 * See the discussion there for explanation of the original problem.
 */
#define JMAX 20
double gsl_integ_simpson(double (*func)(double), double a, double b, double eps)
{
  int j;
  double s, st, ost, os;

  ost = os = -1.e50;

  for(j=1; j<=JMAX; j++){

    trapezoid_rule(func, a, b, j, &st);
    s = (4.*st - ost) / 3.;
    
    if(fabs(s-os) < eps * fabs(os))
      return s;

    os = s;
    ost = st;
  }
  
  GSL_MESSAGE("simpson: too many steps");

  return 0.;
}
#undef JMAX


double gsl_integ_simpson_table(const double * x, const double * y, int n)
{
  int i;
  double result = 0.;
  for(i=0; i<n-1; i++) {
    result += 0.5 * (y[i+1] + y[i]) * (x[i+1] - x[i]);
  }
  return result;
}


/* apparatus for lorentzian variable change to remove pole */
static double (*dummy_f)(double);
static double dummy_x0;
static double dummy_w;
static inline double f_y(double y)
{ return dummy_f(dummy_x0 + dummy_w * tan(dummy_w * y)); }

double gsl_integ_lorenz(double (*f)(double),
                        double x0, double w,
                        double a, double b,
                        double eps)
{
  dummy_f = f;
  dummy_x0 = x0;
  dummy_w = fabs(w);

  if(dummy_w < 10.*sqrt(min_double)) {
    char buff[100];
    sprintf(buff,"lorenz_integ: width w= %g  too small", dummy_w);
    GSL_MESSAGE(buff);
    return 0.;
  }
  else {
    double lower_y = atan((a-x0)/dummy_w) / dummy_w;
    double upper_y = atan((b-x0)/dummy_w) / dummy_w;
    return gsl_integ_simpson(f_y, lower_y, upper_y, eps);
  }
}

#endif /* 0 */
