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
void test_discrete_pdf (double (*f) (void), double (*pdf)(unsigned int), 
			const char *name);

void test_shuffle (void);
void test_choose (void);
double test_beta (void);
double test_beta_pdf (double x);
double test_bernoulli (void);
double test_bernoulli_pdf (unsigned int n);
double test_binomial (void);
double test_binomial_pdf (unsigned int n);
double test_binomial_large (void);
double test_binomial_large_pdf (unsigned int n);
double test_cauchy (void);
double test_cauchy_pdf (double x);
double test_chisq (void);
double test_chisq_pdf (double x);
double test_erlang (void);
double test_erlang_pdf (double x);
double test_exponential (void);
double test_exponential_pdf (double x);
double test_exppow0 (void);
double test_exppow0_pdf (double x);
double test_exppow1 (void);
double test_exppow1_pdf (double x);
double test_exppow1a (void);
double test_exppow1a_pdf (double x);
double test_exppow2 (void);
double test_exppow2_pdf (double x);
double test_exppow2a (void);
double test_exppow2a_pdf (double x);
double test_fdist (void);
double test_fdist_pdf (double x);
double test_flat (void);
double test_flat_pdf (double x);
double test_gamma (void);
double test_gamma_pdf (double x);
double test_gamma1 (void);
double test_gamma1_pdf (double x);
double test_gamma_int (void);
double test_gamma_int_pdf (double x);
double test_gamma_large (void);
double test_gamma_large_pdf (double x);
double test_gaussian (void);
double test_gaussian_pdf (double x);
double test_bivariate_gaussian1 (void);
double test_bivariate_gaussian1_pdf (double x);
double test_bivariate_gaussian2 (void);
double test_bivariate_gaussian2_pdf (double x);
double test_bivariate_gaussian3 (void);
double test_bivariate_gaussian3_pdf (double x);
double test_bivariate_gaussian4 (void);
double test_bivariate_gaussian4_pdf (double x);
double test_gumbel1 (void);
double test_gumbel1_pdf (double x);
double test_gumbel2 (void);
double test_gumbel2_pdf (double x);
double test_ugaussian (void);
double test_ugaussian_pdf (double x);
double test_geometric (void);
double test_geometric_pdf (unsigned int x);
double test_geometric1 (void);
double test_geometric1_pdf (unsigned int x);
double test_hypergeometric1 (void);
double test_hypergeometric1_pdf (unsigned int x);
double test_hypergeometric2 (void);
double test_hypergeometric2_pdf (unsigned int x);
double test_hypergeometric3 (void);
double test_hypergeometric3_pdf (unsigned int x);
double test_hypergeometric4 (void);
double test_hypergeometric4_pdf (unsigned int x);
double test_hypergeometric5 (void);
double test_hypergeometric5_pdf (unsigned int x);
double test_levy1 (void);
double test_levy1_pdf (double x);
double test_levy2 (void);
double test_levy2_pdf (double x);
double test_logistic (void);
double test_logistic_pdf (double x);
double test_lognormal (void);
double test_lognormal_pdf (double x);
double test_negative_binomial (void);
double test_negative_binomial_pdf (unsigned int n);
double test_pascal (void);
double test_pascal_pdf (unsigned int n);
double test_pareto (void);
double test_pareto_pdf (double x);
double test_poisson (void);
double test_poisson_pdf (unsigned int x);
double test_poisson_large (void);
double test_poisson_large_pdf (unsigned int x);
double test_dir2d (void);
double test_dir2d_pdf (double x);
double test_dir3dxy (void);
double test_dir3dxy_pdf (double x);
double test_dir3dyz (void);
double test_dir3dyz_pdf (double x);
double test_dir3dzx (void);
double test_dir3dzx_pdf (double x);
double test_tdist1 (void);
double test_tdist1_pdf (double x);
double test_tdist2 (void);
double test_tdist2_pdf (double x);
double test_laplace (void);
double test_laplace_pdf (double x);
double test_weibull (void);
double test_weibull_pdf (double x);
double test_weibull1 (void);
double test_weibull1_pdf (double x);

gsl_rng *r_global;

int
main (void)
{
  gsl_rng_env_setup() ;
  r_global = gsl_rng_alloc (gsl_rng_default);

#define FUNC(x) x, "gsl_ran_" #x
#define FUNC2(x) x, x ## _pdf, "gsl_ran_" #x
  test_pdf (FUNC2(test_bivariate_gaussian1));
  test_pdf (FUNC2(test_bivariate_gaussian2));
  test_pdf (FUNC2(test_bivariate_gaussian3));
  test_pdf (FUNC2(test_bivariate_gaussian4));

  test_shuffle() ;
  test_choose() ;

  test_moments (FUNC (test_ugaussian), 0.0, 100.0, 0.5);
  test_moments (FUNC (test_ugaussian), -1.0, 1.0, 0.68);
  test_moments (FUNC (test_exponential), 0.0, 1.0, 1- exp(-0.5));
  test_moments (FUNC (test_cauchy), 0.0, 10000.0, 0.5);


  test_pdf (FUNC2(test_beta));
  test_pdf (FUNC2(test_cauchy));
  test_pdf (FUNC2(test_chisq));
  test_pdf (FUNC2(test_erlang));
  test_pdf (FUNC2(test_exponential));

  test_pdf (FUNC2(test_exppow0));
  test_pdf (FUNC2(test_exppow1));
  test_pdf (FUNC2(test_exppow1a));
  test_pdf (FUNC2(test_exppow2));
  test_pdf (FUNC2(test_exppow2a));

  test_pdf (FUNC2(test_fdist));
  test_pdf (FUNC2(test_flat));
  test_pdf (FUNC2(test_gamma));
  test_pdf (FUNC2(test_gamma1));
  test_pdf (FUNC2(test_gamma_int));
  test_pdf (FUNC2(test_gamma_large));
  test_pdf (FUNC2(test_gaussian));
  test_pdf (FUNC2(test_ugaussian));

  test_pdf (FUNC2(test_gumbel1));
  test_pdf (FUNC2(test_gumbel2));
  test_pdf (FUNC2(test_levy1));
  test_pdf (FUNC2(test_levy2));
  test_pdf (FUNC2(test_logistic));
  test_pdf (FUNC2(test_lognormal));
  test_pdf (FUNC2(test_pareto));
  test_pdf (FUNC2(test_tdist1));
  test_pdf (FUNC2(test_tdist2));
  test_pdf (FUNC2(test_laplace));
  test_pdf (FUNC2(test_weibull));
  test_pdf (FUNC2(test_weibull1));

  test_pdf (FUNC2(test_dir2d));
  test_pdf (FUNC2(test_dir3dxy));
  test_pdf (FUNC2(test_dir3dyz));
  test_pdf (FUNC2(test_dir3dzx));

  test_discrete_pdf (FUNC2(test_poisson));
  test_discrete_pdf (FUNC2(test_poisson_large));
  test_discrete_pdf (FUNC2(test_bernoulli));
  test_discrete_pdf (FUNC2(test_binomial));
  test_discrete_pdf (FUNC2(test_binomial_large));
  test_discrete_pdf (FUNC2(test_geometric));
  test_discrete_pdf (FUNC2(test_geometric1));
  test_discrete_pdf (FUNC2(test_hypergeometric1));
  test_discrete_pdf (FUNC2(test_hypergeometric2));
  test_discrete_pdf (FUNC2(test_hypergeometric3));
  test_discrete_pdf (FUNC2(test_hypergeometric4));
  test_discrete_pdf (FUNC2(test_hypergeometric5));
  test_discrete_pdf (FUNC2(test_negative_binomial));
  test_discrete_pdf (FUNC2(test_pascal));

  return gsl_test_summary();
}

void
test_shuffle (void)
{
  double count[10][10] ;
  int x[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9} ;
  int i, j, status = 0;

  for (i = 0; i < 10; i++)
    {
      for (j = 0; j < 10; j++)
	{
	  count[i][j] = 0 ;
	}
    }

  for (i = 0 ; i < N; i++)
    {
      for (j = 0; j < 10; j++)
	x[j] = j ;

      gsl_ran_shuffle (r_global, x, 10, sizeof(int)) ;

      for (j = 0; j < 10; j++)
	count[x[j]][j] ++ ;
    }

  for (i = 0; i < 10; i++)
    {
      for (j = 0; j < 10; j++)
	{
	  double expected = N / 10.0 ;
	  double d = fabs(count[i][j] - expected);
	  double sigma = d / sqrt(expected) ;
	  if (sigma > 5 && d > 1)
	    {
	      status = 1 ;
	      gsl_test (status, 
			"gsl_ran_shuffle %d,%d (%g observed vs %g expected)", 
			i, j, count[i][j]/N, 0.1) ;
	    }
	}
    }
  
  gsl_test (status, "gsl_ran_shuffle on {0, 1, 2, 3, 4, 5, 6, 7, 8, 9}") ;

}

void
test_choose (void)
{
  double count[10] ;
  int x[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9} ;
  int y[3] = {0, 1, 2} ;
  int i, j, status = 0;

  for (i = 0; i < 10; i++)
    {
      count[i] = 0 ;
    }

  for (i = 0 ; i < N; i++)
    {
      for (j = 0; j < 10; j++)
	x[j] = j ;

      gsl_ran_choose (r_global, y, 3, x, 10, sizeof(int)) ;

      for (j = 0; j < 3; j++)
	count[y[j]]++ ;
    }

  for (i = 0; i < 10; i++)
    {
      double expected = 3.0 * N / 10.0 ;
      double d = fabs(count[i] - expected);
      double sigma = d / sqrt(expected) ;
      if (sigma > 5 && d > 1)
	{
	  status = 1 ;
	  gsl_test (status, 
		    "gsl_ran_choose %d (%g observed vs %g expected)", 
		    i, count[i]/N, 0.1) ;
	}
    }
  
  gsl_test (status, "gsl_ran_choose (3) on {0, 1, 2, 3, 4, 5, 6, 7, 8, 9}") ;

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

  gsl_test (status, "%s [%g,%g] (%g observed vs %g expected)",
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
#define STEPS 100
      double sum = 0 ;
      
      if (fabs(x) < 1e-10) /* hit the origin exactly */
	x = 0.0 ; 
      
      for (j = 1; j < STEPS; j++)
	sum += pdf(x + j * dx / STEPS) ;

      p[i] =  0.5 * (pdf(x) + 2*sum + pdf(x + dx - 1e-7)) * dx / STEPS ;
    }

  for (i = 0; i < BINS; i++)
    {
      double x = a + i * dx ;
      double d = fabs(count[i] - N*p[i]) ;
      if (p[i] != 0)
	{
	  double s = d / sqrt(N*p[i]) ;
	  status_i = (s > 5) && (d > 1) ;
	}
      else
	{
	  status_i = (count[i] != 0) ;
	}
      status |= status_i ;
      if (status_i) 
	gsl_test (status_i, "%s [%g,%g) (%g observed vs %g expected)", 
		  name, x, x+dx, count[i]/N, p[i]) ;
    }

  if (status == 0)
    gsl_test (status, "%s, sampling against pdf over range [%g,%g) ", 
	      name, a, b) ;
}

void
test_discrete_pdf (double (*f) (void), double (*pdf)(unsigned int), const char *name)
{
  double count[BINS], p[BINS];
  unsigned int i ;
  int status = 0, status_i =0 ;

  for (i = 0; i < BINS; i++)
    count[i] = 0 ;

  for (i = 0; i < N; i++)
    {
      int r = f ();
      if (r>= 0 && r < BINS)
	count[r]++;
    }
  
  for (i = 0; i < BINS; i++)
    p[i] =  pdf(i) ;

  for (i = 0; i < BINS; i++)
    {
      double d = fabs(count[i] - N*p[i]) ;
      if (p[i] != 0)
	{
	  double s = d/sqrt(N*p[i]) ;
	  status_i = (s > 5) && (d > 1);
	}
      else
	{
	  status_i = (count[i] != 0) ;
	}
      status |= status_i ;
      if (status_i) 
	gsl_test (status_i, "%s i=%d (%g observed vs %g expected)", 
		  name, i, count[i]/N, p[i]) ;
    }

  if (status == 0)
    gsl_test (status, "%s, sampling against pdf over range [%d,%d) ", 
	      name, 0, BINS) ;
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
test_bernoulli (void)
{
  return gsl_ran_bernoulli (r_global, 0.3);
}

double
test_bernoulli_pdf (unsigned int n)
{
  return gsl_ran_bernoulli_pdf (n, 0.3);
}


double
test_binomial (void)
{
  return gsl_ran_binomial (r_global, 0.3, 5);
}

double
test_binomial_pdf (unsigned int n)
{
  return gsl_ran_binomial_pdf (n, 0.3, 5);
}

double
test_binomial_large (void)
{
  return gsl_ran_binomial (r_global, 0.3, 55);
}

double
test_binomial_large_pdf (unsigned int n)
{
  return gsl_ran_binomial_pdf (n, 0.3, 55);
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
  return gsl_ran_chisq (r_global, 13.0);
}

double
test_chisq_pdf (double x)
{
  return gsl_ran_chisq_pdf (x, 13.0);
}

double
test_dir2d (void)
{
  double x=0, y=0, theta;
  gsl_ran_dir_2d (r_global, &x, &y);
  theta = atan2(x,y);
  return theta;
}

double
test_dir2d_pdf (double x)
{
  if (x > -M_PI && x <= M_PI)
    {
      return 1 / (2 * M_PI) ;
    }
  else
    {
      return 0 ;
    }
}


double
test_dir3dxy (void)
{
  double x=0, y=0, z=0, theta;
  gsl_ran_dir_3d (r_global, &x, &y, &z);
  theta = atan2(x,y);
  return theta;
}

double
test_dir3dxy_pdf (double x)
{
  if (x > -M_PI && x <= M_PI)
    {
      return 1 / (2 * M_PI) ;
    }
  else
    {
      return 0 ;
    }
}

double
test_dir3dyz (void)
{
  double x=0, y=0, z=0, theta;
  gsl_ran_dir_3d (r_global, &x, &y, &z);
  theta = atan2(y,z);
  return theta;
}

double
test_dir3dyz_pdf (double x)
{
  if (x > -M_PI && x <= M_PI)
    {
      return 1 / (2 * M_PI) ;
    }
  else
    {
      return 0 ;
    }
}

double
test_dir3dzx (void)
{
  double x=0, y=0, z=0, theta;
  gsl_ran_dir_3d (r_global, &x, &y, &z);
  theta = atan2(z,x);
  return theta;
}

double
test_dir3dzx_pdf (double x)
{
  if (x > -M_PI && x <= M_PI)
    {
      return 1 / (2 * M_PI) ;
    }
  else
    {
      return 0 ;
    }
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
test_exppow0 (void)
{
  return gsl_ran_exppow (r_global, 3.7, 0.3);
}

double
test_exppow0_pdf (double x)
{
  return gsl_ran_exppow_pdf (x, 3.7, 0.3);
}

double
test_exppow1 (void)
{
  return gsl_ran_exppow (r_global, 3.7, 1.0);
}

double
test_exppow1_pdf (double x)
{
  return gsl_ran_exppow_pdf (x, 3.7, 1.0);
}

double
test_exppow1a (void)
{
  return gsl_ran_exppow (r_global, 3.7, 1.9);
}

double
test_exppow1a_pdf (double x)
{
  return gsl_ran_exppow_pdf (x, 3.7, 1.9);
}

double
test_exppow2 (void)
{
  return gsl_ran_exppow (r_global, 3.7, 2.0);
}

double
test_exppow2_pdf (double x)
{
  return gsl_ran_exppow_pdf (x, 3.7, 2.0);
}


double
test_exppow2a (void)
{
  return gsl_ran_exppow (r_global, 3.7, 7.5);
}

double
test_exppow2a_pdf (double x)
{
  return gsl_ran_exppow_pdf (x, 3.7, 7.5);
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
test_gamma1 (void)
{
  return gsl_ran_gamma (r_global, 1.0);
}

double
test_gamma1_pdf (double x)
{
  return gsl_ran_gamma_pdf (x, 1.0);
}


double
test_gamma_int (void)
{
  return gsl_ran_gamma (r_global, 10.0);
}

double
test_gamma_int_pdf (double x)
{
  return gsl_ran_gamma_pdf (x, 10.0);
}


double
test_gamma_large (void)
{
  return gsl_ran_gamma (r_global, 20.0);
}

double
test_gamma_large_pdf (double x)
{
  return gsl_ran_gamma_pdf (x, 20.0);
}


double
test_gaussian (void)
{
  return gsl_ran_gaussian (r_global, 3.0);
}

double
test_gaussian_pdf (double x)
{
  return gsl_ran_gaussian_pdf (x, 3.0);
}

double
test_bivariate_gaussian1 (void)
{
  double x = 0, y = 0;
  gsl_ran_bivariate_gaussian (r_global, 3.0, 2.0, 0.3, &x, &y);
  return x ;
}

double
test_bivariate_gaussian1_pdf (double x)
{
  return gsl_ran_gaussian_pdf (x, 3.0);
}

double
test_bivariate_gaussian2 (void)
{
  double x = 0, y = 0;
  gsl_ran_bivariate_gaussian (r_global, 3.0, 2.0, 0.3, &x, &y);
  return y ;
}

double
test_bivariate_gaussian2_pdf (double y)
{
  int i, n = 10 ;
  double sum = 0 ;
  double a = -10, b = 10, dx = (b - a)/n ;
  for (i = 0; i < n ; i++)
    {
      double x = a + i * dx ;
      sum += gsl_ran_bivariate_gaussian_pdf (x, y, 3.0, 2.0, 0.3) * dx ;
    }
  return sum ;
}


double
test_bivariate_gaussian3 (void)
{
  double x = 0, y = 0;
  gsl_ran_bivariate_gaussian (r_global, 3.0, 2.0, 0.3, &x, &y);
  return x + y ;
}

double
test_bivariate_gaussian3_pdf (double x)
{
  double sx = 3.0, sy = 2.0, r = 0.3;
  double su = (sx+r*sy) ;
  double sv = sy*sqrt(1-r*r) ;
  double sigma = sqrt(su*su + sv*sv) ;
    
  return gsl_ran_gaussian_pdf (x, sigma);
}

double
test_bivariate_gaussian4 (void)
{
  double x = 0, y = 0;
  gsl_ran_bivariate_gaussian (r_global, 3.0, 2.0, -0.5, &x, &y);
  return x + y ;
}

double
test_bivariate_gaussian4_pdf (double x)
{
  double sx = 3.0, sy = 2.0, r = -0.5;
  double su = (sx+r*sy) ;
  double sv = sy*sqrt(1-r*r) ;
  double sigma = sqrt(su*su + sv*sv) ;
    
  return gsl_ran_gaussian_pdf (x, sigma);
}




double
test_ugaussian (void)
{
  return gsl_ran_ugaussian (r_global);
}

double
test_ugaussian_pdf (double x)
{
  return gsl_ran_ugaussian_pdf (x);
}


double
test_geometric (void)
{
  return gsl_ran_geometric (r_global, 0.5);
}

double
test_geometric_pdf (unsigned int n)
{
  return gsl_ran_geometric_pdf (n, 0.5);
}

double
test_geometric1 (void)
{
  return gsl_ran_geometric (r_global, 1.0);
}

double
test_geometric1_pdf (unsigned int n)
{
  return gsl_ran_geometric_pdf (n, 1.0);
}

double
test_hypergeometric1 (void)
{
  return gsl_ran_hypergeometric (r_global, 5, 7, 4);
}

double
test_hypergeometric1_pdf (unsigned int n)
{
  return gsl_ran_hypergeometric_pdf (n, 5, 7, 4);
}


double
test_hypergeometric2 (void)
{
  return gsl_ran_hypergeometric (r_global, 5, 7, 11);
}

double
test_hypergeometric2_pdf (unsigned int n)
{
  return gsl_ran_hypergeometric_pdf (n, 5, 7, 11);
}

double
test_hypergeometric3 (void)
{
  return gsl_ran_hypergeometric (r_global, 5, 7, 1);
}

double
test_hypergeometric3_pdf (unsigned int n)
{
  return gsl_ran_hypergeometric_pdf (n, 5, 7, 1);
}

double
test_hypergeometric4 (void)
{
  return gsl_ran_hypergeometric (r_global, 5, 7, 20);
}

double
test_hypergeometric4_pdf (unsigned int n)
{
  return gsl_ran_hypergeometric_pdf (n, 5, 7, 20);
}

double
test_hypergeometric5 (void)
{
  return gsl_ran_hypergeometric (r_global, 2, 7, 5);
}

double
test_hypergeometric5_pdf (unsigned int n)
{
  return gsl_ran_hypergeometric_pdf (n, 2, 7, 5);
}



double
test_gumbel1 (void)
{
  return gsl_ran_gumbel1 (r_global, 3.12, 4.56);
}

double
test_gumbel1_pdf (double x)
{
  return gsl_ran_gumbel1_pdf (x, 3.12, 4.56);
}

double
test_gumbel2 (void)
{
  return gsl_ran_gumbel2 (r_global, 3.12, 4.56);
}

double
test_gumbel2_pdf (double x)
{
  return gsl_ran_gumbel2_pdf (x, 3.12, 4.56);
}

double
test_levy1 (void)
{
  return gsl_ran_levy (r_global, 5.0, 1.0);
}

double
test_levy1_pdf (double x)
{
  return gsl_ran_cauchy_pdf (x, 5.0);
}

double
test_levy2 (void)
{
  return gsl_ran_levy (r_global, 5.0, 2.0);
}

double
test_levy2_pdf (double x)
{
  return gsl_ran_gaussian_pdf (x, sqrt(2.0) * 5.0 );
}


double
test_logistic (void)
{
  return gsl_ran_logistic (r_global, 3.1);
}

double
test_logistic_pdf (double x)
{
  return gsl_ran_logistic_pdf (x, 3.1);
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
test_negative_binomial (void)
{
  return gsl_ran_negative_binomial (r_global, 0.3, 20.0);
}

double
test_negative_binomial_pdf (unsigned int n)
{
  return gsl_ran_negative_binomial_pdf (n, 0.3, 20.0);
}

double
test_pascal (void)
{
  return gsl_ran_pascal (r_global, 0.8, 3);
}

double
test_pascal_pdf (unsigned int n)
{
  return gsl_ran_pascal_pdf (n, 0.8, 3);
}


double
test_pareto (void)
{
  return gsl_ran_pareto (r_global, 1.9, 2.75);
}

double
test_pareto_pdf (double x)
{
  return gsl_ran_pareto_pdf (x, 1.9, 2.75);
}

double
test_poisson (void)
{
  return gsl_ran_poisson (r_global, 5.0);
}

double
test_poisson_pdf (unsigned int n)
{
  return gsl_ran_poisson_pdf (n, 5.0);
}

double
test_poisson_large (void)
{
  return gsl_ran_poisson (r_global, 30.0);
}

double
test_poisson_large_pdf (unsigned int n)
{
  return gsl_ran_poisson_pdf (n, 30.0);
}


double
test_tdist1 (void)
{
  return gsl_ran_tdist (r_global, 1.75);
}

double
test_tdist1_pdf (double x)
{
  return gsl_ran_tdist_pdf (x, 1.75);
}

double
test_tdist2 (void)
{
  return gsl_ran_tdist (r_global, 12.75);
}

double
test_tdist2_pdf (double x)
{
  return gsl_ran_tdist_pdf (x, 12.75);
}


double
test_laplace (void)
{
  return gsl_ran_laplace (r_global, 2.75);
}

double
test_laplace_pdf (double x)
{
  return gsl_ran_laplace_pdf (x, 2.75);
}

double
test_weibull (void)
{
  return gsl_ran_weibull (r_global, 3.14, 2.75);
}

double
test_weibull_pdf (double x)
{
  return gsl_ran_weibull_pdf (x, 3.14, 2.75);
}


double
test_weibull1 (void)
{
  return gsl_ran_weibull (r_global, 2.97, 1.0);
}

double
test_weibull1_pdf (double x)
{
  return gsl_ran_weibull_pdf (x, 2.97, 1.0);
}
