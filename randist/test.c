
#include <config.h>
#include <stdio.h>
#include <stdlib.h>

#include <gsl_ran.h>
#include <gsl_randist.h>
#include <gsl_test.h>

double
test_poisson ()
{
  return gsl_ran_poisson (10);
};

test_rndist (double (*f) (), double mean, double var, int n)
{
  int i;
  double sum = 0;
  double est_mean = 0, est_var = 0, est_skew = 0, est_four = 0;

  for (i = 0; i < n; ++i)
    {
      const double r = f ();
      const double delta = r - mean;
      est_mean += (r - est_mean) / (i + 1);
      est_var += (delta * delta - est_var) / (i + 1);
      est_skew += (delta * delta * delta - est_skew) / (i + 1);
      est_four += (delta * delta * delta * delta - est_four) / (i + 1);
    }
  printf ("est_mean = %g\n", est_mean);
  printf ("est_var = %g\n", est_var);
  printf ("est_skew = %g\n", est_skew);
  printf ("est_four = %g\n", est_four);
};


int
main (int argc, char **argv)
{
  int i, n = 100000;
  int randseed = 17;

  gsl_ran_seed (randseed);

  test_rndist (gsl_ran_gaussian, 0, 1, n);
  test_rndist (test_poisson, 10, sqrt (10), n);


  return 0;
}
