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

double beta (void);
double beta_pdf (double x);
double cauchy (void);
double cauchy_pdf (double x);
double chisq (void);
double chisq_pdf (double x);
double erlang (void);
double erlang_pdf (double x);
double exponential (void);
double exponential_pdf (double x);
double fdist (void);
double fdist_pdf (double x);
double flat (void);
double flat_pdf (double x);
double gamma (void);
double gamma_pdf (double x);
double gaussian (void);
double gaussian_pdf (double x);
double geometric (void);
double logistic (void);
double logistic_pdf (double x);
double lognormal (void);
double lognormal_pdf (double x);
double pareto (void);
double pareto_pdf (double x);
double poisson (void);
double tdist (void);
double tdist_pdf (double x);
double weibull (void);
double weibull_pdf (double x);

gsl_rng *r_global;

int
main (void)
{
  r_global = gsl_rng_alloc (gsl_rng_default);

#define FUNC(x) x, "gsl_ran_" #x
  test_moments (FUNC (gaussian), 0.0, 100.0, 0.5);
  test_moments (FUNC (gaussian), -1.0, 1.0, 0.68);
  test_moments (FUNC (exponential), 0.0, 1.0, 1- exp(-0.5));
  test_moments (FUNC (cauchy), 0.0, 10000.0, 0.5);

#define FUNC2(x) x, x ## _pdf, "gsl_ran_" #x
  test_pdf (FUNC2(beta));
  test_pdf (FUNC2(cauchy));
  test_pdf (FUNC2(chisq));
  test_pdf (FUNC2(erlang));
  test_pdf (FUNC2(exponential));
  test_pdf (FUNC2(fdist));
  test_pdf (FUNC2(flat));
  test_pdf (FUNC2(gamma));
  test_pdf (FUNC2(gaussian));
  test_pdf (FUNC2(logistic));
  test_pdf (FUNC2(lognormal));
  test_pdf (FUNC2(pareto));
  test_pdf (FUNC2(tdist));
  test_pdf (FUNC2(weibull));

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
beta (void)
{
  return gsl_ran_beta (r_global, 2.0, 3.0);
}

double
beta_pdf (double x)
{
  return gsl_ran_beta_pdf (x, 2.0, 3.0);
}

double
cauchy (void)
{
  return gsl_ran_cauchy (r_global, 2.0);
}

double
cauchy_pdf (double x)
{
  return gsl_ran_cauchy_pdf (x, 2.0);
}

double
chisq (void)
{
  return gsl_ran_chisq (r_global, 3.0);
}

double
chisq_pdf (double x)
{
  return gsl_ran_chisq_pdf (x, 3.0);
}


double
erlang (void)
{
  return gsl_ran_erlang (r_global, 3.0, 4.0);
}

double
erlang_pdf (double x)
{
  return gsl_ran_erlang_pdf (x, 3.0, 4.0);
}

double
exponential (void)
{
  return gsl_ran_exponential (r_global, 2.0);
}

double
exponential_pdf (double x)
{
  return gsl_ran_exponential_pdf (x, 2.0);
}

double
fdist (void)
{
  return gsl_ran_fdist (r_global, 3.0, 4.0);
}

double
fdist_pdf (double x)
{
  return gsl_ran_fdist_pdf (x, 3.0, 4.0);
}

double
flat (void)
{
  return gsl_ran_flat (r_global, 3.0, 4.0);
}

double
flat_pdf (double x)
{
  return gsl_ran_flat_pdf (x, 3.0, 4.0);
}

double
gamma (void)
{
  return gsl_ran_gamma (r_global, 2.5);
}

double
gamma_pdf (double x)
{
  return gsl_ran_gamma_pdf (x, 2.5);
}

double
gaussian (void)
{
  return gsl_ran_gaussian (r_global);
}

double
gaussian_pdf (double x)
{
  return gsl_ran_gaussian_pdf (x);
}


double
geometric (void)
{
  return gsl_ran_geometric (r_global, 0.5);
}

double
geometric_pdf (unsigned int x)
{
  return gsl_ran_geometric_pdf (x, 0.5);
}


double
logistic (void)
{
  return gsl_ran_logistic (r_global);
}

double
logistic_pdf (double x)
{
  return gsl_ran_logistic_pdf (x);
}

double
lognormal (void)
{
  return gsl_ran_lognormal (r_global);
}

double
lognormal_pdf (double x)
{
  return gsl_ran_lognormal_pdf (x);
}


double
pareto (void)
{
  return gsl_ran_pareto (r_global, 2.75);
}

double
pareto_pdf (double x)
{
  return gsl_ran_pareto_pdf (x, 2.75);
}

double
poisson (void)
{
  return gsl_ran_poisson (r_global, 2.0);
}

double
tdist (void)
{
  return gsl_ran_tdist (r_global, 2.75);
}

double
tdist_pdf (double x)
{
  return gsl_ran_tdist_pdf (x, 2.75);
}

double
weibull (void)
{
  return gsl_ran_weibull (r_global, 2.75);
}

double
weibull_pdf (double x)
{
  return gsl_ran_weibull_pdf (x, 2.75);
}
