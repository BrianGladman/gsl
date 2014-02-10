/* multifit/test_nonlinear.c
 * 
 * Copyright (C) 2007, 2013, 2014 Brian Gough, Patrick Alken
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

typedef struct
{
  const char *name;
  double *x0;       /* initial parameters (size p) */
  double *x_sol;    /* solution parameters (size p) */
  double *f_sumsq;  /* ||f(x_sol)||^2 */
  double *sigma;
  double *epsrel;   /* relative tolerance for solution checking */
  gsl_multifit_function_fdf *fdf;
} test_fdf_problem;

#include "test_brown.c"
#include "test_fn.c"

#include "test_bard.c"
#include "test_enso.c"
#include "test_hahn1.c"
#include "test_helical.c"
#include "test_kirby2.c"
#include "test_lin1.c"
#include "test_powell1.c"
#include "test_powell2.c"
#include "test_rosenbrock.c"
#include "test_roth.c"

static void test_lmder (gsl_multifit_function_fdf * f, double x0[], 
                        double * X, double F[], double * cov);
static void test_fdf(const gsl_multifit_fdfsolver_type * T, const double xtol,
                     const double gtol, const double ftol, test_fdf_problem *problem);

/*
 * These test problems are taken from
 *
 * H. B. Nielsen, UCTP test problems for unconstrained optimization,
 * IMM Department of Mathematical Modeling, Tech. Report IMM-REP-2000-17,
 * 2000.
 */
static test_fdf_problem *test_fdf_nielsen[] = {
  &lin1_problem,
  &rosenbrock_problem,
  &helical_problem,
  &powell1_problem,
  &powell2_problem,
  &roth_problem,
  &bard_problem,

  NULL
};

/* NIST test cases */
static test_fdf_problem *test_fdf_nist[] = {
  &enso_problem,
  &kirby2_problem,
  &hahn1_problem,
  NULL
};

static void
test_nonlinear(void)
{
  const double xtol = 1e-15;
  const double gtol = 1e-15;
  const double ftol = 0.0;
  gsl_multifit_function_fdf f;
  size_t i;

  for (i = 0; test_fdf_nielsen[i] != NULL; ++i)
    {
      test_fdf(gsl_multifit_fdfsolver_lmniel, xtol, gtol, ftol,
               test_fdf_nielsen[i]);
    }

  for (i = 0; test_fdf_nist[i] != NULL; ++i)
    {
      test_fdf(gsl_multifit_fdfsolver_lmniel, xtol, gtol, ftol,
               test_fdf_nist[i]);
      test_fdf(gsl_multifit_fdfsolver_lmsder, 1e-5, 1e-5, 0.0, test_fdf_nist[i]);
      test_fdf(gsl_multifit_fdfsolver_lmder, 1e-5, 1e-5, 0.0, test_fdf_nist[i]);
    }

  exit(1);

  {
    f = make_fdf (&brown_f, &brown_df, &brown_fdf, brown_N, brown_P, 0);
    test_lmder(&f, brown_x0, &brown_X[0][0], brown_F, &brown_cov[0][0]);
  }
}

static void
test_lmder (gsl_multifit_function_fdf * f, double x0[], 
            double * X, double F[], double * cov)
{
  const gsl_multifit_fdfsolver_type *T;
  gsl_multifit_fdfsolver *s;

  const size_t n = f->n;
  const size_t p = f->p;

  size_t iter = 0, i;
  
  gsl_vector_view x = gsl_vector_view_array (x0, p);

  T = gsl_multifit_fdfsolver_lmsder;
  s = gsl_multifit_fdfsolver_alloc (T, n, p);
  gsl_multifit_fdfsolver_set (s, f, &x.vector);

  do
    {
      gsl_multifit_fdfsolver_iterate (s);

      for (i = 0 ; i < p; i++)
        {
          gsl_test_rel (gsl_vector_get (s->x, i), X[p*iter+i], 1e-5, 
                        "lmsder, iter=%u, x%u", iter, i);
        }

      gsl_test_rel (gsl_blas_dnrm2 (s->f), F[iter], 1e-5, 
                    "lmsder, iter=%u, f", iter);

      iter++;
    }
  while (iter < 20);
  
  {
    size_t i, j;
    gsl_matrix * covar = gsl_matrix_alloc (4, 4);
    gsl_multifit_covar (s->J, 0.0, covar);

    for (i = 0; i < 4; i++) 
      {
        for (j = 0; j < 4; j++)
          {
            gsl_test_rel (gsl_matrix_get(covar,i,j), cov[i*p + j], 1e-7, 
                          "gsl_multifit_covar cov(%d,%d)", i, j) ;
          }
      }

    gsl_matrix_free (covar);
  }

  gsl_multifit_fdfsolver_free (s);

}

static void
test_fdf(const gsl_multifit_fdfsolver_type * T, const double xtol,
         const double gtol, const double ftol, test_fdf_problem *problem)
{
  gsl_multifit_function_fdf *fdf = problem->fdf;
  const size_t n = fdf->n;
  const size_t p = fdf->p;
  const double *sigma = problem->sigma;
  const double epsrel = *(problem->epsrel);
  const size_t max_iter = 500;
  gsl_vector_view x0 = gsl_vector_view_array(problem->x0, p);
  gsl_multifit_fdfsolver *s = gsl_multifit_fdfsolver_alloc (T, n, p);
  const char *pname = problem->name;
  const char *sname = gsl_multifit_fdfsolver_name(s);
  size_t iter = 0;
  int status, info;
  size_t i;

  gsl_multifit_fdfsolver_set(s, fdf, &x0.vector);

  printf("working on %s/%s\n", sname, pname);

  do
    {
      status = gsl_multifit_fdfsolver_iterate (s);
      gsl_test(status, "%s/%s iterate status=%s", sname, pname,
               gsl_strerror(status));

      status = gsl_multifit_test_convergence(s, xtol, gtol, ftol, &info);
    }
  while (status == GSL_CONTINUE && ++iter < max_iter);

  gsl_test(status, "%s/%s did not converge", sname, pname);

  printf("iter = %zu\n", iter);

  /* check computed x = x_sol */
  for (i = 0; i < p; ++i)
    {
      double xi = gsl_vector_get(s->x, i);
      double xi_exact = problem->x_sol[i];

      gsl_test_rel(xi, xi_exact, epsrel,
                   "%s/%s solution i=%zu",
                   sname, pname, i);
    }

  {
    double s2;
    double s2_exact = *(problem->f_sumsq);

    /* check computed/exact ||f||^2 */
    gsl_blas_ddot(s->f, s->f, &s2);
    gsl_test_rel(s2, s2_exact, epsrel,
                 "%s/%s sumsq", sname, pname);

    /* check variances */
    if (sigma)
      {
        gsl_matrix * covar = gsl_matrix_alloc (p, p);
        gsl_multifit_covar (s->J, 0.0, covar);

        for (i = 0; i < p; i++) 
          {
            double ei = sqrt(s2/(n-p))*sqrt(gsl_matrix_get(covar,i,i));
            gsl_test_rel (ei, sigma[i], epsrel, 
                          "%s/%s, sigma(%d)", sname, pname, i) ;
          }

        gsl_matrix_free (covar);
      }
  }

  gsl_multifit_fdfsolver_free(s);
}
