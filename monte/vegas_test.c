/* Little routine to test vegas integration. */

/* Author: M.J. Booth */
/* RCS: $Id$ */

/* Use the examples in Lepage's paper: */

#include <config.h>
#include <math.h>
#include <stdio.h>

#include <gsl_errno.h>
#include <gsl_test.h>

#include <gsl_vegas.h>

#define SQR(x)  ((x)*(x))

double a = 0.1;
double c = 0;
int num_dim = 1;

double xl[11] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
double xu[11] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

double f0(double x[]);
double f1(double x[]);
double f2(double x[]);
double f3(double x[]);

int main() 
{
  double res = 0;
  double err = 0;
  double chisq = 0;
  int status = 0;
  double tol = 1e-2;
  int step = 1;

  c = (1.0 + sqrt(10.0))/9.0 ;

  mode = 1;
  /*  verbose = 1; */
  calls = 1000;
  max_it_num = 10;

  gsl_test_verbose(1);

  printf("Testing product function\n");
  for (num_dim = 1; num_dim < 10; num_dim++) {
    gsl_monte_vegas(f0, xl, xu, num_dim, &res, &err, &chisq);
    gsl_test_rel(res, 1.0, tol, "vegas(f0), dim=%d, err=%.4f, chisq=%.4f", 
		 num_dim, err, chisq); 
  }

  mode = 1;
  /*  verbose = 1;*/
  calls = 10000;
  alpha = 2.0;
  max_it_num = 10;
  tol = 1e-2;

  printf("Testing single gaussian\n");
  for (num_dim = 1; num_dim < 10; num_dim++) {
    if (num_dim > 4) 
      alpha -= 0.2; /* Lepage has=1 for n=9 so we step toward that. */
    if (num_dim > 8)
      calls = 15000;
    status = gsl_monte_vegas(f1, xl, xu, num_dim, &res, &err, &chisq);
    gsl_test_rel(res, 1.0, tol, "vegas(f1), dim=%d, err=%.4f, chisq=%.4f", 
		 num_dim, err, chisq); 
  }

  mode = 1;
  /*  verbose = 1; */
  calls = 20000;
  alpha = 2.0;
  max_it_num = 5;
  tol = 2e-2;
  
  step = 2;
  printf("Testing double gaussian\n");
  for (num_dim = 1; num_dim < 10; num_dim += step) {
    switch (num_dim) {
    case 5:
      calls = 40000;
      break;
    case 6:
      calls = 50000;
      break;
    case 7:
      max_it_num = 10;
      tol = 0.08;
      calls = 70000;
      break;
    case 9:
      calls = 100000;
      max_it_num = 20;
      tol = 0.04;
      alpha = 1.0;
      break;
    }
    /* let's try same trick of stepping alpha */
    if (num_dim > 4) 
      alpha -= (double) step/10.0; 
    status = gsl_monte_vegas(f2, xl, xu, num_dim, &res, &err, &chisq);
    gsl_test_rel(res, 1.0, tol, "vegas(f2), dim=%d, err=%.4f, chisq=%.4f", 
		 num_dim, err, chisq); 
  }

  mode = 1;
  /*  verbose = 1; */
  calls = 10000;
  alpha = 1.5;
  max_it_num = 10;
  tol = 1e-2;

  printf("Testing Tsuda's function\n");
  for (num_dim = 1; num_dim < 10; num_dim++) {
    status = gsl_monte_vegas(f3, xl, xu, num_dim, &res, &err, &chisq);
    gsl_test_rel(res, 1.0, tol, "vegas(f3), dim=%d, err=%.4f, chisq=%.4f", 
		 num_dim, err, chisq); 
  }

  return gsl_test_summary();
}

/* Simple product integral */
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
  
