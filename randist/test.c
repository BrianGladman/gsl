#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl_randist.h>
#include <gsl_rng.h>
#include <gsl_test.h>

#define N 100000
void test_moments (double (*f) (void), const char *name,
		   double a, double b, double p);
void test_pdf (double (*f) (void), double (*pdf)(double), const char *name);

double test_beta (void);
double test_beta_pdf (double x);
double test_cauchy (void);
double test_cauchy_pdf (double x);
double test_chisq (void);
double test_chisq_pdf (double x);
double test_erlang (void);
double test_erlang_pdf (double x);
double test_exponential (void);
double test_exponential_pdf (double x);
double test_fdist (void);
double test_fdist_pdf (double x);
double test_flat (void);
double test_flat_pdf (double x);
double test_gamma (void);
double test_gamma_pdf (double x);
double test_gaussian (void);
double test_gaussian_pdf (double x);
double test_geometric (void);
double test_logistic (void);
double test_logistic_pdf (double x);
double test_lognormal (void);
double test_lognormal_pdf (double x);
double test_pareto (void);
double test_pareto_pdf (double x);
double test_poisson (void);
double test_tdist (void);
double test_tdist_pdf (double x);
double test_weibull (void);
double test_weibull_pdf (double x);

gsl_rng *r_global;

int
main (void)
{
  r_global = gsl_rng_alloc (gsl_rng_default);

#define FUNC(x) x, "gsl_ran_" #x
  test_moments (FUNC (test_gaussian), 0.0, 100.0, 0.5);
  test_moments (FUNC (test_gaussian), -1.0, 1.0, 0.68);
  test_moments (FUNC (test_exponential), 0.0, 1.0, 1- exp(-0.5));
  test_moments (FUNC (test_cauchy), 0.0, 10000.0, 0.5);

#define FUNC2(x) x, x ## _pdf, "gsl_ran_" #x
  test_pdf (FUNC2(test_beta));
  test_pdf (FUNC2(test_cauchy));
  test_pdf (FUNC2(test_chisq));
  test_pdf (FUNC2(test_erlang));
  test_pdf (FUNC2(test_exponential));
  test_pdf (FUNC2(test_fdist));
  test_pdf (FUNC2(test_flat));
  test_pdf (FUNC2(test_gamma));
  test_pdf (FUNC2(test_gaussian));
  test_pdf (FUNC2(test_logistic));
  test_pdf (FUNC2(test_lognormal));
  test_pdf (FUNC2(test_pareto));
  test_pdf (FUNC2(test_tdist));
  test_pdf (FUNC2(test_weibull));

  return 0;
}

void
test_moments (double (*f) (void), const char *name,
	      double a, double b, double p)
{
  int i;
  double count = 0, expected, sigma;
  int status;

  for (i = 0; i < N; i++)
    {
      double r = f ();
      if (r < b && r > a)
	count++;
    }

  expected = p * N;
  sigma = fabs (count - expected) / sqrt (expected);

  status = (sigma > 3);

  gsl_test (status, "%s, range [%g,%g] (%g observed vs %g expected)",
	    name, a, b, count / N, p);
}

#define BINS 100

void
test_pdf (double (*f) (void), double (*pdf)(double), const char *name)
{
  double count[BINS], p[BINS];
  double a = -5.0, b = +5.0 ;
  double dx = (b - a) / BINS ;
  int i,j,status = 0, status_i =0 ;

  for (i = 0; i < BINS; i++)
    count[i] = 0 ;

  for (i = 0; i < N; i++)
    {
      double r = f ();
      if (r < b && r > a)
	{ 
	  j =  (r - a)/dx ;
	  count[j]++;
	}
    }
  
  for (i = 0; i < BINS; i++)
    {
      /* Compute an approximation to the integral of p(x) from x to
         x+dx using Simpson's rule */

      double x = a + i * dx ;

      double sum = 0 ;
      for (j = 1; j < 1000; j++)
	sum += pdf(x + j * dx / 1000) ;

      p[i] =  0.5 * (pdf(x) + 2*sum + pdf(x + dx - 1e-7)) * dx / 1000 ;
    }

  for (i = 0; i < BINS; i++)
    {
      double x = a + i * dx ;
      double d = fabs(count[i] - N*p[i]) ;
      if (p[i] != 0)
	{
	  d = d / sqrt(N*p[i]) ;
	  status_i = (d > 5) ;
	}
      else
	{
	  status_i = (count[i] != 0) ;
	}
      status |= status_i ;
      if (status_i) 
	gsl_test (status_i, "%s, range [%g,%g) (%g observed vs %g expected)", 
		  name, x, x+dx, count[i]/N, p[i]) ;
    }

  if (status == 0)
    gsl_test (status, "%s, sampling against pdf, range [%g,%g) ", name, a, b) ;
}
      

double
test_beta (void)
{
  return gsl_ran_beta (r_global, 2.0, 3.0);
}

double
test_beta_pdf (double x)
{
  return gsl_ran_beta_pdf (x, 2.0, 3.0);
}

double
test_cauchy (void)
{
  return gsl_ran_cauchy (r_global, 2.0);
}

double
test_cauchy_pdf (double x)
{
  return gsl_ran_cauchy_pdf (x, 2.0);
}

double
test_chisq (void)
{
  return gsl_ran_chisq (r_global, 3.0);
}

double
test_chisq_pdf (double x)
{
  return gsl_ran_chisq_pdf (x, 3.0);
}


double
test_erlang (void)
{
  return gsl_ran_erlang (r_global, 3.0, 4.0);
}

double
test_erlang_pdf (double x)
{
  return gsl_ran_erlang_pdf (x, 3.0, 4.0);
}

double
test_exponential (void)
{
  return gsl_ran_exponential (r_global, 2.0);
}

double
test_exponential_pdf (double x)
{
  return gsl_ran_exponential_pdf (x, 2.0);
}

double
test_fdist (void)
{
  return gsl_ran_fdist (r_global, 3.0, 4.0);
}

double
test_fdist_pdf (double x)
{
  return gsl_ran_fdist_pdf (x, 3.0, 4.0);
}

double
test_flat (void)
{
  return gsl_ran_flat (r_global, 3.0, 4.0);
}

double
test_flat_pdf (double x)
{
  return gsl_ran_flat_pdf (x, 3.0, 4.0);
}

double
test_gamma (void)
{
  return gsl_ran_gamma (r_global, 2.5);
}

double
test_gamma_pdf (double x)
{
  return gsl_ran_gamma_pdf (x, 2.5);
}

double
test_gaussian (void)
{
  return gsl_ran_gaussian (r_global);
}

double
test_gaussian_pdf (double x)
{
  return gsl_ran_gaussian_pdf (x);
}


double
test_geometric (void)
{
  return gsl_ran_geometric (r_global, 0.5);
}

double
test_geometric_pdf (unsigned int x)
{
  return gsl_ran_geometric_pdf (x, 0.5);
}

double
test_logistic (void)
{
  return gsl_ran_logistic (r_global);
}

double
test_logistic_pdf (double x)
{
  return gsl_ran_logistic_pdf (x);
}

double
test_lognormal (void)
{
  return gsl_ran_lognormal (r_global);
}

double
test_lognormal_pdf (double x)
{
  return gsl_ran_lognormal_pdf (x);
}

double
test_pareto (void)
{
  return gsl_ran_pareto (r_global, 2.75);
}

double
test_pareto_pdf (double x)
{
  return gsl_ran_pareto_pdf (x, 2.75);
}

double
test_poisson (void)
{
  return gsl_ran_poisson (r_global, 2.0);
}

double
test_tdist (void)
{
  return gsl_ran_tdist (r_global, 2.75);
}

double
test_tdist_pdf (double x)
{
  return gsl_ran_tdist_pdf (x, 2.75);
}

double
test_weibull (void)
{
  return gsl_ran_weibull (r_global, 2.75);
}

double
test_weibull_pdf (double x)
{
  return gsl_ran_weibull_pdf (x, 2.75);
}
