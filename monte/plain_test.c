/* monte/plain_test.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Michael Booth
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/* Little routine to test plain MC integration. */

/* Author: M.J. Booth */
/* RCS: $Id$ */

/* Use the examples in Lepage's paper: */

#include <config.h>
#include <math.h>
#include <stdio.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_monte_plain.h>

#define SQR(x)  ((x)*(x))

double a = 0.1;
double c = 0;
int num_dim = 1;

double xl[11] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
double xu[11] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
double xu2[11] = {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};
double xu3[2] = {GSL_DBL_MAX, GSL_DBL_MAX};

double fconst(double x[]);
double f0(double x[]);
double f1(double x[]);
double f2(double x[]);
double f3(double x[]);

void my_error_handler (const char *reason, const char *file,
                       int line, int err);

int main() 
{
  double res = 0;
  double err = 0;
  /* double chisq = 0; */
  int status = 0;
  double tol = 1e-2;
  int step = 1;
  int calls = 1000;

  gsl_monte_plain_state* s = gsl_monte_plain_alloc(10);

  c = (1.0 + sqrt(10.0))/9.0 ;

  calls = 1000;

  gsl_test_verbose(1);

  gsl_set_error_handler (&my_error_handler);

  gsl_ieee_env_setup (); 

  printf("testing allocation/input checks\n");

  status = gsl_monte_plain_validate(s, xl, xu, 4, 10);
  gsl_test(status != 0,  "error if not initialized");

  gsl_monte_plain_init(s);
  status = gsl_monte_plain_validate(s, xl, xu3, 1, 1);
  gsl_test(status != 0,  "error if limits too large");
  status = gsl_monte_plain_validate(s, xl, xu, 0, 10);
  gsl_test(status != 0,  "error if num_dim = 0");
  status = gsl_monte_plain_validate(s, xl, xu, 1, 0);
  gsl_test(status != 0,  "error if calls = 0");
  status = gsl_monte_plain_validate(s, xu, xl, 1, 10);
  gsl_test(status != 0,  "error if xu < xl");


  printf("Testing constant function and normalization\n");
  for (num_dim = 1; num_dim < 10; num_dim++) {
    calls *= 1.2;
    status = gsl_monte_plain_integrate(s, fconst, xl, xu2, num_dim, calls, 
				       &res, &err);
    gsl_test_rel(res, pow(2, num_dim), tol, 
		 "plain(fconst), calls=%d, dim=%d, err=%.4f", 
		 calls, num_dim, err); 
  }

  calls = 3000;
  tol = 2e-2;
  printf("Testing product function\n");
  for (num_dim = 1; num_dim < 10; num_dim++) {
    calls *= 1.8;
    status = gsl_monte_plain_integrate(s, f0, xl, xu, num_dim, calls, 
				       &res, &err);
    gsl_test_rel(res, 1.0, tol, "plain(f0), calls=%d, dim=%d, err=%.4f", 
		 calls, num_dim, err); 
  }

  /*  verbose = 1;*/
  calls =  20000;
  tol = 0.04;
  step = 2;

  printf("Testing single gaussian\n");
  for (num_dim = 1; num_dim < 8; num_dim += step) {
    switch (num_dim) {
    case 3:
      calls = 40000;
      tol = 0.1;
      break;
    case 5:
      calls = 120000;
      tol = 0.15;
      break;
    case 6:
      calls = 250000;
      tol = 0.1;
      break;
    case 7:
      calls = 500000;
      tol = 0.15;
      break;
    case 8:
      tol = 0.1;
      break;
    case 9:
      calls = 1000000;
      tol = 0.14;
      break;
    }
    status = gsl_monte_plain_integrate(s, f1, xl, xu, num_dim, calls, 
				       &res, &err);
    gsl_test_rel(res, 1.0, tol, "plain(f1), calls=%d, dim=%d, err=%.4f", 
		 calls, num_dim, err); 
  }

  /*  verbose = 1; */
  calls = 20000;
  tol = 2e-2;
  step = 2;

  printf("Testing double gaussian\n");
  for (num_dim = 1; num_dim < 10; num_dim += step) {
    switch (num_dim) {
    case 3:
      tol = 0.09;
      calls = 40000;
      break;
    case 5:
      tol = 0.15;
      calls = 100000;
      break;
    case 7:
      calls = 120000;
      tol = 0.15;
      break;
    case 9:
      calls = 1000000;
      tol = 0.5; /* lazy... */
      break;
    }
    status = gsl_monte_plain_integrate(s, f2, xl, xu, num_dim, calls, 
				       &res, &err);
    gsl_test_rel(res, 1.0, tol, "plain(f2), calls=%d, dim=%d, err=%.4f", 
		 calls, num_dim, err); 
  }

  /*  verbose = 1; */
  calls = 10000;
  tol = 4e-2;

  printf("Testing Tsuda's function\n");
  for (num_dim = 1; num_dim < 10; num_dim++) {
    calls *= 1.2;
    if ( num_dim == 6 )
      tol *= 3;
    status = gsl_monte_plain_integrate(s, f3, xl, xu, num_dim, calls, 
				       &res, &err);
    gsl_test_rel(res, 1.0, tol, "plain(f3), calls=%d, dim=%d, err=%.4f", 
		 calls, num_dim, err); 
  }

  return gsl_test_summary();
}

/* Simple constant function */
double fconst(double x[])
{
  return  1;
}

/* Simple product function */
double f0(double x[])
{
  double prod = 1.0;
  int i;

  for (i=0; i < num_dim; ++i) {
    prod *= 2.0*x[i];
  }

  return  prod;
}

/* Gaussian centered at 1/2. */

double f1(double x[])
{
  double sum = 0.;

  int i;
  for (i=0; i < num_dim; i++) {
    sum += SQR(x[i] - 0.5);
  }
  return ( pow(M_2_SQRTPI/(2.*a), (double)num_dim) * exp(- sum/SQR(a)) );
}

/* double gaussian */
double f2(double x[])
{
  double sum1 = 0.;
  double sum2 = 0.;

  int i;
  for (i=0; i < num_dim; i++) {
    sum1 += SQR(x[i] - 1./3.);
    sum2 += SQR(x[i] - 2./3.);
  }
  return 0.5*pow(M_2_SQRTPI/(2.*a), num_dim)*( exp(- sum1/SQR(a)) + exp(-sum2/SQR(a)) );
}



/* Tsuda's example */
double f3(double x[])
{
  double prod = 1.;

  int i;
  for (i=0; i < num_dim; i++) {
    prod *= c/(c+1)*SQR((c+1)/(c+x[i]));
  }

  return prod;
}
  

void
my_error_handler (const char *reason, const char *file, int line, int err)
{
  if (0) printf ("(caught [%s:%d: %s (%d)])\n", file, line, reason, err) ;
}
