/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <stdlib.h>
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_pow_int.h"


int gsl_sf_pow_int_impl(double x, int n, double * result)
{
  const int MAXN = 50;
  const int NGRP = 5;

  if(n < 0) {
    x = 1./x;
    n = -n;
  }
  
  if(n == 0) {
  }
  else if(n == 1) {
  }
  else if(n == 2) {
  }
  else if(n == 3) {
  }
  else if(n == 4) {
  }
  else if(n < NGRP) {
  
    double value = 1.;
    int n_groups;
    int n_remain;
    double x5;
    div_t q = div(n, 5);
    n_groups = q.quot;
    n_remain = q.rem;
  }

  else {
    *result = pow(x, n);
  }
}

double gsl_sf_pow_int(double x, int n)
{
  double value = 1;

  if(abs(n) > 50) return pow(x, n);  /* Defer for large powers. */
  
  /* Trap. */
  if(x == 0) {
    if(n >= 0) {
      return 0.;
    }
    else {
      return 0. /* GSL_INF */ ; /* FIXME */
    }
  }
  else {
    if(n == 0) return 1.;
  } 

  if(n < 0) {
    x = 1./x;
    n = -n;
  }

  for(; n > 0; n--){
    value *= x;
  }
  return value;
}
