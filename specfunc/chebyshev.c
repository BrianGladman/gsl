#include <stdlib.h>
#include <math.h>
#include "constants.h"
#include "memory.h"
#include "error.h"
#include "chebyshev.h"


struct ChebSeries * cheb_calc(double (*func)(double),
			      double a, double b,
			      int order)
{
  int k, j;
  double bma = 0.5 * (b - a);
  double bpa = 0.5 * (b + a);
  double fac = 2./(order +1.);
  double * f = new_vector_d(order+1);

  struct ChebSeries * cs = (struct ChebSeries *)
    malloc(sizeof(struct ChebSeries));
  if(cs == 0) {
    char buff[100];
    sprintf(buff, "cheb_calc: allocation failure");
    push_error(buff, Error_Alloc_);
    return 0;
  }

  cs->order = order;
  cs->a = a;
  cs->b = b;
  cs->c = new_vector_d(order+1);

  for(k = 0; k<=order; k++) {
    double y = cos(constPi_ * (k+0.5)/(order+1));
    f[k] = func(y*bma + bpa);
  }

  for(j = 0; j<=order; j++) {
    double sum = 0.0;
    for(k = 0; k<=order; k++) sum += f[k]*cos(constPi_ * j*(k+0.5)/(order+1));
    cs->c[j] = fac * sum;
  }

  free_vector(f);

  return cs;
}


double cheb_eval(double x, const struct ChebSeries * cs)
{
  int j;
  double d  = 0.;
  double dd = 0.;

  double y  = (2.*x - cs->a - cs->b) / (cs->b - cs->a);
  double y2 = 2. * y;

  for(j = cs->order; j>=1; j--) {
    double temp = d;
    d = y2*d - dd + cs->c[j];
    dd = temp;
  }
  return y*d - dd + 0.5 * cs->c[0];
}


void free_cheb(struct ChebSeries * cs)
{
  if(cs != 0) {
    if(cs->c != 0) free_vector(cs);
    free(cs);
  }
}
