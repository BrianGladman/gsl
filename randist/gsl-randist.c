#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <gsl_randist.h>
#include <gsl_rng.h>
#include <gsl_test.h>

void error (const char * s);


int
main (int argc, char *argv[])
{
  size_t i;
  size_t n = 0;
  double mu = 0, nu = 0, nu1 = 0, nu2 = 0, sigma = 0, a = 0, b = 0, p = 0;
  double sigmax = 0, sigmay = 0, rho = 0;
  double x = 0, y =0, z=0  ;
  unsigned int N = 0, t = 0, n1 = 0, n2 = 0 ;
  unsigned long int seed = 0 ;
  const char * name ;
  gsl_rng * r ;

  if (argc < 4) 
    {
      printf (
"Usage: gsl-randist seed n DIST param1 param2 ...\n"
"Generates n samples from the distribution DIST with parameters param1,\n"
"param2, etc. Valid distributions are,\n"
"\n"
"  exponential\n"
"  levy\n") ;
      exit (0);
    }

  argv++ ; seed = atol (argv[0]); argc-- ;
  argv++ ; n = atol (argv[0]); argc-- ;
  argv++ ; name = argv[0] ; argc-- ; argc-- ;

  gsl_rng_env_setup() ;

  if (gsl_rng_default_seed != 0) {
    fprintf(stderr, 
	    "overriding GSL_RNG_SEED with command line value, seed = %ld\n", 
	    seed) ;
  }
  
  gsl_rng_default_seed = seed ;

  r = gsl_rng_alloc(gsl_rng_default) ;


#define NAME(x) !strcmp(name,(x))
#define OUTPUT(x) for (i = 0; i < n; i++) { printf("%g\n", (x)) ; }
#define OUTPUT2(a,x,y) for(i = 0; i < n; i++) { a ; printf("%g %g\n", x, y) ; }
#define OUTPUT3(a,x,y,z) for(i = 0; i < n; i++) { a ; printf("%g %g %g\n", x, y, z) ; }
#define INT_OUTPUT(x) for (i = 0; i < n; i++) { printf("%d\n", (x)) ; }
#define ARGS(x,y) if (argc != x) error(y) ;
#define DBL_ARG(x) if (argc) { x=atof((++argv)[0]);argc--;} else {error( #x);};
#define INT_ARG(x) if (argc) { x=atoi((++argv)[0]);argc--;} else {error( #x);};

  if (NAME("bernoulli"))
    {
      ARGS(1, "p = probability of success");
      DBL_ARG(p)
      INT_OUTPUT(gsl_ran_bernoulli (r, p));
    }
  else if (NAME("beta"))
    {
      ARGS(2, "a,b = shape parameters");
      DBL_ARG(a)
      DBL_ARG(b)
      OUTPUT(gsl_ran_beta (r, a, b));
    }
  else if (NAME("binomial"))
    {
      ARGS(2, "p = probability, N = number of trials");
      DBL_ARG(p)
      INT_ARG(N)
      INT_OUTPUT(gsl_ran_binomial (r, p, N));
    }
  else if (NAME("cauchy"))
    {
      ARGS(1, "mu = scale parameter");
      DBL_ARG(mu)
      OUTPUT(gsl_ran_cauchy (r, mu));
    }
  else if (NAME("chisq"))
    {
      ARGS(1, "nu = degrees of freedom");
      DBL_ARG(nu)
      OUTPUT(gsl_ran_chisq (r, nu));
    }
  else if (NAME("erlang"))
    {
      ARGS(2, "a = scale parameter, b = order");
      DBL_ARG(a)
      DBL_ARG(b)
      OUTPUT(gsl_ran_erlang (r, a, b));
    }
  else if (NAME("exponential"))
    {
      ARGS(1, "mu = mean value");
      DBL_ARG(mu) ;
      OUTPUT(gsl_ran_exponential (r, mu));
    }
  else if (NAME("exppow"))
    {
      ARGS(2, "mu = scale parameter, a = power (1=exponential, 2=gaussian)");
      DBL_ARG(mu) ;
      DBL_ARG(a) ;
      OUTPUT(gsl_ran_exppow (r, mu, a));
    }
  else if (NAME("fdist"))
    {
      ARGS(2, "nu1, nu2 = degrees of freedom parameters");
      DBL_ARG(nu1) ;
      DBL_ARG(nu2) ;
      OUTPUT(gsl_ran_fdist (r, nu1, nu2));
    }
  else if (NAME("flat"))
    {
      ARGS(2, "a = lower limit, b = upper limit");
      DBL_ARG(a) ;
      DBL_ARG(b) ;
      OUTPUT(gsl_ran_flat (r, a, b));
    }
  else if (NAME("gamma"))
    {
      ARGS(2, "a = order, b = scale");
      DBL_ARG(a) ;
      DBL_ARG(b) ;
      OUTPUT(gsl_ran_gamma (r, a, b));
    }
  else if (NAME("gaussian"))
    {
      ARGS(1, "sigma = standard deviation");
      DBL_ARG(sigma) ;
      OUTPUT(gsl_ran_gaussian (r, sigma));
    }
  else if (NAME("bivariate-gaussian"))
    {
      ARGS(3, "sigmax = x std.dev., sigmay = y std.dev., rho = correlation");
      DBL_ARG(sigmax) ;
      DBL_ARG(sigmay) ;
      DBL_ARG(rho) ;
      OUTPUT2(gsl_ran_bivariate_gaussian (r, sigmax, sigmay, rho, &x, &y), 
	      x, y);
    }
  else if (NAME("dir-2d"))
    {
      OUTPUT2(gsl_ran_dir_2d (r, &x, &y), x, y);
    }
  else if (NAME("dir-3d"))
    {
      OUTPUT3(gsl_ran_dir_3d (r, &x, &y, &z), x, y, z);
    }
  else if (NAME("geometric"))
    {
      ARGS(1, "p = bernoulli trial probability of success");
      DBL_ARG(p) ;
      INT_OUTPUT(gsl_ran_geometric (r, p));
    }
  else if (NAME("gumbel1"))
    {
      ARGS(2, "a = order, b = scale parameter");
      DBL_ARG(a) ;
      DBL_ARG(b) ;
      OUTPUT(gsl_ran_gumbel1 (r, a, b));
    }
  else if (NAME("gumbel2"))
    {
      ARGS(2, "a = order, b = scale parameter");
      DBL_ARG(a) ;
      DBL_ARG(b) ;
      OUTPUT(gsl_ran_gumbel2 (r, a, b));
    }
  else if (NAME("hypergeometric"))
    {
      ARGS(3, "n1 = tagged population, n2 = untagged population, t = number of trials");
      INT_ARG(n1) ;
      INT_ARG(n2) ;
      INT_ARG(t) ;
      INT_OUTPUT(gsl_ran_hypergeometric (r, n1, n2, t));
    }
  else if (NAME("laplace"))
    {
      ARGS(1, "mu = scale parameter");
      DBL_ARG(mu) ;
      OUTPUT(gsl_ran_laplace (r, mu));
    }
  else if (NAME("levy"))
    {
      ARGS(2, "mu = scale parameter, a = power (1=cauchy, 2=gaussian)");
      DBL_ARG(mu) ;
      DBL_ARG(a) ;
      OUTPUT(gsl_ran_levy (r, mu, a));
    }
  else if (NAME("logarithmic"))
    {
      ARGS(1, "p = probability");
      DBL_ARG(p) ;
      INT_OUTPUT(gsl_ran_logarithmic (r, p));
    }
  else if (NAME("logistic"))
    {
      ARGS(1, "mu = scale parameter");
      DBL_ARG(mu) ;
      OUTPUT(gsl_ran_logistic (r, mu));
    }
  else if (NAME("lognormal"))
    {
      ARGS(2, "zeta = location parameter, sigma = scale parameter");
      DBL_ARG(mu) ;
      DBL_ARG(sigma) ;
      OUTPUT(gsl_ran_exppow (r, mu, sigma));
    }
  else if (NAME("negative-binomial"))
    {
      ARGS(2, "p = probability, a = order");
      DBL_ARG(p) ;
      DBL_ARG(a) ;
      INT_OUTPUT(gsl_ran_negative_binomial (r, p, a));
    }
  else if (NAME("pareto"))
    {
      ARGS(2, "mu = scale parameter, a = power");
      DBL_ARG(mu) ;
      DBL_ARG(a) ;
      OUTPUT(gsl_ran_pareto (r, mu, a));
    }
  else if (NAME("pascal"))
    {
      ARGS(2, "p = probability, n = order (integer)");
      DBL_ARG(mu) ;
      INT_ARG(N) ;
      INT_OUTPUT(gsl_ran_pascal (r, mu, N));
    }
  else if (NAME("poisson"))
    {
      ARGS(1, "mu = scale parameter");
      DBL_ARG(mu) ;
      INT_OUTPUT(gsl_ran_poisson (r, mu));
    }
  else if (NAME("tdist"))
    {
      ARGS(1, "nu = degrees of freedom");
      DBL_ARG(nu) ;
      OUTPUT(gsl_ran_tdist (r, nu));
    }
  else if (NAME("weibull"))
    {
      ARGS(2, "mu = scale parameter, a = exponent");
      DBL_ARG(mu) ;
      DBL_ARG(a) ;
      OUTPUT(gsl_ran_weibull (r, mu, a));
    }
  else
    {
      fprintf(stderr,"Error: unrecognized distribution: %s\n", name) ;
    }

  return 0 ;
}


void
error (const char * s)
{
  fprintf(stderr, "Error: arguments should be %s\n",s) ;
  exit (EXIT_FAILURE) ;
}
