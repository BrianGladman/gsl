#include <config.h>
#include <stdio.h>
#include <stdlib.h>

#include <gsl_randist.h>
#include <gsl_rng.h>
#include <gsl_test.h>

#define N 10000
void test_moments (double (*f) (void), double mean, double var);
double gaussian (void);
double exponential (void);

gsl_rng * r_global ;

int
main (void)
{
  r_global = gsl_rng_alloc (gsl_rng_default) ;

  test_moments (gaussian, 0.0, 1.0);
  test_moments (exponential, 0.0, 1.0);

  return 0;
}

void
test_moments (double (*f) (void), double mean, double var)
{
  int i;
  double est_mean = 0, est_var = 0, est_skew = 0, est_four = 0;

  for (i = 0; i < N; i++)
    {
      const double r = f ();
      const double delta = r - mean;
      est_mean += (r - est_mean) / (i + 1);
      est_var += (delta * delta - est_var) / (i + 1);
      est_skew += (delta * delta * delta - est_skew) / (i + 1);
      est_four += (delta * delta * delta * delta - est_four) / (i + 1);
    }
  printf ("est_mean = %g vs %g\n", est_mean, mean);
  printf ("est_var = %g vs %g\n", est_var, var);
  printf ("est_skew = %g\n", est_skew);
  printf ("est_four = %g\n", est_four);
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

