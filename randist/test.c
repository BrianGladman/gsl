#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl_randist.h>
#include <gsl_rng.h>
#include <gsl_test.h>

#define N 100000
void test_moments (double (*f) (void), const char * name, 
		   double a, double b, double p);

double beta (void);
double cauchy (void);
double chisq (void);
double erlang (void);
double exponential (void);
double fdist (void);
double flat (void);
double gamma (void);
double gaussian (void);
double geometric (void);
double logistic (void);
double lognormal (void);
double pareto (void);
double poisson (void);
double tdist (void);
double weibull (void);

gsl_rng * r_global ;

int
main (void)
{
  r_global = gsl_rng_alloc (gsl_rng_default) ;

#define FUNC(x) x, "gsl_ran_" #x
  test_moments (FUNC(gaussian), 0.0, 100.0, 0.5);
  test_moments (FUNC(gaussian), -1.0, 1.0, 0.68);

  test_moments (FUNC(exponential), 0.0, 1.0, 0.63212);
  test_moments (FUNC(cauchy), 0.0, 10000.0, 0.5);

  return 0;
}

void
test_moments (double (*f) (void), const char * name, 
	      double a, double b, double p)
{
  int i ;
  double count = 0, expected, sigma ;
  int status ;

  for (i = 0; i < N; i++)
    {
      double r = f ();
      if (r < b && r > a)
	count++ ;
    }
  
  expected = p * N ;
  sigma = fabs(count - expected) / sqrt(expected) ;
  
  status = (sigma > 3) ;

  gsl_test(status, "%s, range [%g,%g] (%g observed vs %g expected)",
	  name, a, b, count/N, p) ;
}

double 
beta (void)
{
  return gsl_ran_beta (r_global, 1.0, 2.0) ;
}

double 
cauchy (void)
{
  return gsl_ran_cauchy (r_global, 1.0) ;
}

double 
chisq (void)
{
  return gsl_ran_chisq (r_global, 1.0) ;
}

double 
erlang (void)
{
  return gsl_ran_erlang (r_global, 1.0, 2.0) ;
}

double 
exponential (void)
{
  return gsl_ran_exponential (r_global, 1.0) ;
}

double 
fdist (void)
{
  return gsl_ran_fdist (r_global, 1.0, 2.0) ;
}

double 
flat (void)
{
  return gsl_ran_flat (r_global, 1.0, 2.0) ;
}

double 
gamma (void)
{
  return gsl_ran_gamma (r_global, 1.0) ;
}

double 
gaussian (void)
{
  return gsl_ran_gaussian (r_global) ;
}

double 
geometric (void)
{
  return gsl_ran_geometric (r_global, 0.5) ;
}

double 
logistic (void)
{
  return gsl_ran_logistic (r_global) ;
}

double 
lognormal (void)
{
  return gsl_ran_lognormal (r_global) ;
}

double 
pareto (void)
{
  return gsl_ran_pareto (r_global, 2.0) ;
}

double 
poisson (void)
{
  return gsl_ran_poisson (r_global, 2.0) ;
}

double 
tdist (void)
{
  return gsl_ran_tdist (r_global, 2.0) ;
}

double 
weibull (void)
{
  return gsl_ran_weibull (r_global, 2.0) ;
}



