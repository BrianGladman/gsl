#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_integration.h>

#include "testint.c"

gsl_integration_workspace *w;
int inverse;


double
cauchy (double x, void *p)
{
  return gsl_ran_cauchy_pdf (x, *(double*)p);
}

double
gaussian (double x, void *p)
{
  return gsl_ran_gaussian_pdf (x, *(double*)p);
}

double
laplace (double x, void *p)
{
  return gsl_ran_laplace_pdf (x, *(double*)p);
}

double
rayleigh (double x, void *p)
{
  return gsl_ran_rayleigh_pdf (x, *(double*)p);
}

double
flat (double x, void *p)
{
  double * c = (double *)p;
  return gsl_ran_flat_pdf (x, c[0], c[1]);
}

double
lognormal (double x, void *p)
{
  double * c = (double *)p;
  return gsl_ran_lognormal_pdf (x, c[0], c[1]);
}

double
gamma (double x, void *p)
{
  double * c = (double *)p;
  return gsl_ran_gamma_pdf (x, c[0], c[1]);
}

double
chisq (double x, void *p)
{
  double * c = (double *)p;
  return gsl_ran_chisq_pdf (x, c[0]);
}

double
fdist (double x, void *p)
{
  double * c = (double *)p;
  return gsl_ran_fdist_pdf (x, c[0], c[1]);
}

double
tdist (double x, void *p)
{
  double * c = (double *)p;
  return gsl_ran_tdist_pdf (x, c[0]);
}

double
beta (double x, void *p)
{
  double * c = (double *)p;
  return gsl_ran_beta_pdf (x, c[0], c[1]);
}

double
gumbel1 (double x, void *p)
{
  double * c = (double *)p;
  return gsl_ran_gumbel1_pdf (x, c[0], c[1]);
}

double
gumbel2 (double x, void *p)
{
  double * c = (double *)p;
  return gsl_ran_gumbel2_pdf (x, c[0], c[1]);
}

double
weibull (double x, void *p)
{
  double * c = (double *)p;
  return gsl_ran_weibull_pdf (x, c[0], c[1]);
}

double
pareto (double x, void *p)
{
  double * c = (double *)p;
  return gsl_ran_pareto_pdf (x, c[0], c[1]);
}

double
logistic (double x, void *p)
{
  double * c = (double *)p;
  return gsl_ran_logistic_pdf (x, c[0]);
}

double
integrate (gsl_function * f, double a, double b)
{
  double res = 0, err = 0;
  double tol = 1e-12;
  gsl_integration_qag (f, a, b, 0, tol, 1000, 0,  w, &res, &err);
  return res;
}

double
integrate_lower (gsl_function * f, double b)
{
  double res = 0, err = 0;
  qagil (f, b, 0, 1e-12, 1000, w, &res, &err);
  return res;
}

double
integrate_upper (gsl_function * f, double a)
{
  double res = 0, err = 0;
  qagiu (f, a, 0, 1e-12, 1000, w, &res, &err);
  return res;
}


int
test (const char * name, gsl_function * f, double x[], int N, char *fmt)
{
  int i;
  double res, err, sum = 0, sumerr = 0;

  printf ("void test_auto_%s (void);\n\n", name);

  printf ("void\ntest_auto_%s (void)\n{\n", name);

  /* gsl_set_error_handler_off(); */

  w = gsl_integration_workspace_alloc (1000);

  for (i = 0; i < N; i++)
    {
      res = 0;
      err = 0;

      if (x[0] < -1000)
        {
          if (x[i] < 0)
            {
              res = integrate_lower (f, x[i]);
            }
          else
            {
              res = 1 - integrate_upper (f, x[i]);
            }
          sum = res;
        }
      else
        {
          if (i == 0)
            sum += 0;
          else
            sum += integrate(f, x[i-1], x[i]);
        }
      
      if (res < 0) 
        continue;

      printf (fmt, "_P", x[i], sum);

      if (inverse && (sum != 0 && sum != 1) && (x[i] == 0 || sum * 1e-4 < GSL_FN_EVAL(f,x[i]) * fabs(x[i])))
        printf (fmt, "_Pinv", sum, x[i]);
    }

  printf("\n");

  sum=0;
  sumerr=0;

  for (i = N-1; i >= 0; i--)
    {
      res = 0;
      err = 0;

      if (x[N-1] > 1000)
        {
          if (x[i] > 0)
            {
              res = integrate_upper (f, x[i]);
            }
          else
            {
              res = 1-integrate_lower (f, x[i]);
            }

          sum = res;
        }
      else 
        {
          if (i == N-1)
            sum += 0;
          else
            sum += integrate(f, x[i], x[i+1]);
        }

      printf (fmt, "_Q", x[i], sum);

      if (inverse && (sum != 0 && sum != 1) && (x[i] == 0 || sum * 1e-4 < GSL_FN_EVAL(f,x[i]) * fabs(x[i])))
        printf (fmt, "_Qinv", sum, x[i]);
    }
  
  printf ("}\n\n");

  gsl_integration_workspace_free (w);

}

int
main (void)
{
  const int N = 43;

  double xall[] = {
    -1e10, -1e9, -1e8, -1e7, -1e6, -1e5, -1e4, -1e3, -1e2, -1e1,
    -1,
    -1e-1, -1e-2, -1e-3, -1e-4, -1e-5, -1e-6, -1e-7, -1e-8, -1e-9, -1e-10,
    0,
    1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1,
    1,
    1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10
  };

  const int Npos = 22;

  double xpos[] = {
    0,
    1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1,
    1,
    1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10
  };

  const int N01 = 23;

  double x01[] = {
    0,
    1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1,
    0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 0.9,
    0.99, 0.999, 0.9999, 0.99999, 
    1.0
  };

#define TEST(name,params,range,n) { double p[] = { params } ; gsl_function f = {&name, p}; test(#name, &f, range, n, "    TEST(gsl_cdf_" #name "%s, (%.16e," #params "), %.12e, TEST_TOL6);\n"); }

#define TEST2(name,p1,p2,range,n) { double p[] = { p1,p2 } ; gsl_function f = {&name, p}; test(#name, &f, range, n, "    TEST(gsl_cdf_" #name "%s, (%.16e," #p1 "," #p2 "), %.12e, TEST_TOL6);\n"); }

#define TEST2A(desc,name,p1,p2,range,n) { double p[] = { p1,p2 } ; gsl_function f = {&name, p}; test(#desc, &f, range, n, "    TEST(gsl_cdf_" #name "%s, (%.16e," #p1 "," #p2 "), %.12e, TEST_TOL6);\n"); }


  inverse = 0;
  TEST2(beta, 1.3, 2.7, x01, N01);
  TEST2(fdist, 5.3, 2.7, xpos, Npos);

  inverse = 1;
  TEST(cauchy, 1.3, xall, N);
  TEST(gaussian, 1.3, xall, N);
  TEST(laplace, 1.3, xall, N);
  TEST(rayleigh, 1.3, xpos, Npos);

  TEST2(flat, 1.3, 750.0, xpos, Npos);
  TEST2(lognormal, 1.3,2.7, xpos, Npos);
  TEST2(gamma, 1.3,2.7, xpos, Npos);

  TEST(chisq, 1.3, xpos, Npos);
  TEST(tdist, 1.3, xall, N);

  TEST2(gumbel1, 1.3, 2.7, xall, N);
  TEST2(gumbel2, 1.3, 2.7, xpos, Npos);
  TEST2(weibull, 1.3, 2.7, xpos, Npos);

  TEST2(pareto, 1.3, 2.7, xpos, Npos);
  TEST(logistic, 1.3, xall, N);

  TEST2A(gammalarge, gamma, 1.3,123.0, xpos, Npos);
}
