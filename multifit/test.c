/* multifit/test.c
 * 
 * Copyright (C) 2007, 2013 Brian Gough, Patrick Alken
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

/* These tests are based on the NIST Statistical Reference Datasets
   See http://www.nist.gov/itl/div898/strd/index.html for more
   information. */

#include <config.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>

#include <gsl/gsl_ieee_utils.h>

#include "test_longley.c"
#include "test_filip.c"
#include "test_pontius.c"
#include "test_brown.c"
#include "test_enso.c"
#include "test_kirby2.c"
#include "test_hahn1.c"
#include "test_nelson.c"
#include "test_fn.c"
#include "test_estimator.c"

void
test_lmder (gsl_multifit_function_fdf * f, double x0[], 
            double * X, double F[], double * cov);

void
test_fdf (const char * name, gsl_multifit_function_fdf * f, 
          double x0[], double x[], double sumsq,
          double sigma[], double epsrel, double epsrel_sigma);

void test_ridge(void);

int
main (void)
{
  double epsrel = 1.0e-5;
  double epsrel_sigma = 1.0e-4;
  double epsrel_fd = 1.0e-2;
  double epsrel_sigma_fd = 1.0e-2;
  gsl_multifit_function_fdf f;

  gsl_ieee_env_setup();

  /* test linear regression */
  test_longley();
  test_filip();
  test_pontius();
  test_estimator();
  test_ridge();

  {
    f = make_fdf (&brown_f, &brown_df, &brown_fdf, brown_N, brown_P, 0);
    test_lmder(&f, brown_x0, &brown_X[0][0], brown_F, &brown_cov[0][0]);
  }

  {
    f = make_fdf (&enso_f, &enso_df, &enso_fdf, enso_N, enso_P, 0);
    test_fdf("nist-ENSO", &f, enso_x0, enso_x, enso_sumsq, enso_sigma, epsrel, epsrel_sigma);

    f = make_fdf (&enso_f, NULL, NULL, enso_N, enso_P, 0);
    test_fdf("nist-ENSO", &f, enso_x0, enso_x, enso_sumsq, enso_sigma, epsrel_fd, epsrel_sigma_fd);

  }

  {
    f = make_fdf (&kirby2_f, &kirby2_df, &kirby2_fdf, kirby2_N, kirby2_P, 0);
    test_fdf("nist-kirby2", &f, kirby2_x0, kirby2_x, kirby2_sumsq, kirby2_sigma, epsrel, epsrel_sigma);

    f = make_fdf (&kirby2_f, NULL, NULL, kirby2_N, kirby2_P, 0);
    test_fdf("nist-kirby2", &f, kirby2_x0, kirby2_x, kirby2_sumsq, kirby2_sigma, epsrel_fd, epsrel_sigma_fd);
  }

  {
    f = make_fdf (&hahn1_f, &hahn1_df, &hahn1_fdf, hahn1_N, hahn1_P, 0);
    test_fdf("nist-hahn1", &f, hahn1_x0, hahn1_x, hahn1_sumsq, hahn1_sigma, epsrel, epsrel_sigma);

    f = make_fdf (&hahn1_f, NULL, NULL, hahn1_N, hahn1_P, 0);
    test_fdf("nist-hahn1", &f, hahn1_x0, hahn1_x, hahn1_sumsq, hahn1_sigma, epsrel_fd, epsrel_sigma_fd);
  }

#ifdef JUNK
  {
    f = make_fdf (&nelson_f, &nelson_df, &nelson_fdf, nelson_N, nelson_P, 0);
    test_fdf("nist-nelson", &f, nelson_x0, nelson_x, nelson_sumsq, nelson_sigma, epsrel, epsrel_sigma);
  }
#endif

  /* now summarize the results */

  exit (gsl_test_summary ());
}


void
test_lmder (gsl_multifit_function_fdf * f, double x0[], 
            double * X, double F[], double * cov)
{
  const gsl_multifit_fdfsolver_type *T;
  gsl_multifit_fdfsolver *s;

  const size_t n = f->n;
  const size_t p = f->p;

  int status;
  size_t iter = 0, i;
  
  gsl_vector_view x = gsl_vector_view_array (x0, p);

  T = gsl_multifit_fdfsolver_lmsder;
  s = gsl_multifit_fdfsolver_alloc (T, n, p);
  gsl_multifit_fdfsolver_set (s, f, &x.vector);

  do
    {
      status = gsl_multifit_fdfsolver_iterate (s);

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

void
test_fdf (const char * name, gsl_multifit_function_fdf * f, 
          double x0[], double x_final[], 
          double f_sumsq, double sigma[], double epsrel,
          double epsrel_sigma)
{
  const gsl_multifit_fdfsolver_type *T;
  gsl_multifit_fdfsolver *s;
  
  const size_t n = f->n;
  const size_t p = f->p;

  int status;
  size_t iter = 0;

  gsl_vector_view x = gsl_vector_view_array (x0, p);

  T = gsl_multifit_fdfsolver_lmsder;
  s = gsl_multifit_fdfsolver_alloc (T, n, p);
  gsl_multifit_fdfsolver_set (s, f, &x.vector);

  do
    {
      status = gsl_multifit_fdfsolver_iterate (s);

#ifdef DEBUG
       printf("iter = %d  status = %d  |f| = %.18e x = \n", 
         iter, status, gsl_blas_dnrm2 (s->f));
         
         gsl_vector_fprintf(stdout, s->x, "%.8e");
#endif       
      status = gsl_multifit_test_delta (s->dx, s->x, 0.0, 1e-7);

      iter++;
    }
  while (status == GSL_CONTINUE && iter < 1000);
  
  {
    size_t i;
    gsl_matrix * covar = gsl_matrix_alloc (p, p);
    gsl_multifit_covar (s->J, 0.0, covar);

    for (i = 0 ; i < p; i++)
      {
        gsl_test_rel (gsl_vector_get (s->x, i), x_final[i], epsrel, 
                      "%s, lmsder, x%u", name, i);
      }


    {
      double s2 = pow(gsl_blas_dnrm2 (s->f), 2.0);

      gsl_test_rel (s2, f_sumsq, epsrel, "%s, lmsder, |f|^2", name);

      for (i = 0; i < p; i++) 
        {
          double ei = sqrt(s2/(n-p))*sqrt(gsl_matrix_get(covar,i,i));
          gsl_test_rel (ei, sigma[i], epsrel_sigma, 
                        "%s, sigma(%d)", name, i) ;
        }
    }

    gsl_matrix_free (covar);
  }

  /* check higher level driver routine */
  {
    size_t i;

    gsl_multifit_fdfsolver_set (s, f, &x.vector);
    gsl_multifit_fdfsolver_driver (s, 1000, 0.0, 1.0e-7);

    for (i = 0 ; i < p; i++)
      {
        gsl_test_rel (gsl_vector_get (s->x, i), x_final[i], epsrel, 
                      "%s, lmsder, x%u", name, i);
      }
  }

  /* Check that there is no hidden state, restarting should 
     produce identical results. */

  {
    int status0, status1;
    size_t i;
    gsl_multifit_fdfsolver *t = gsl_multifit_fdfsolver_alloc (T, n, p);
    gsl_multifit_fdfsolver_set (t, f, &x.vector);

    /* do a few extra iterations to stir things up */

    gsl_multifit_fdfsolver_set (s, f, &x.vector);

    for (i = 0; i < 3; i++) 
      {
        gsl_multifit_fdfsolver_iterate (s);
      }

    gsl_multifit_fdfsolver_set (s, f, &x.vector);

    do
      {
        status0 = gsl_multifit_fdfsolver_iterate (s);
        status1 = gsl_multifit_fdfsolver_iterate (t);

        gsl_test_int(status0, status1, "%s, lmsder status after set iter=%u", name, iter);
        
        for (i = 0; i < p; i++) {
          double sxi = gsl_vector_get(s->x,i);
          double txi = gsl_vector_get(t->x,i);
#ifdef DEBUG
          printf("%d %g %g\n", i, sxi, txi);
#endif
          gsl_test_rel(sxi, txi, 1e-15, "%s, lmsder after set, %u/%u", name, iter, i);
        }
        
#ifdef DEBUG
        printf("iter = %d  status = %d  |f| = %.18e x = \n", 
               iter, status, gsl_blas_dnrm2 (s->f));
        
        gsl_vector_fprintf(stdout, s->x, "%.8e");
#endif       
        status0 = gsl_multifit_test_delta (s->dx, s->x, 0.0, 1e-7);
        status1 = gsl_multifit_test_delta (t->dx, s->x, 0.0, 1e-7);
        
        gsl_test_int(status0, status1, "%s, lmsder test delta status after set iter=%u", name, iter);

        iter++;
      }
    while (status1 == GSL_CONTINUE && iter < 1000);

    gsl_multifit_fdfsolver_free (t);
  }

  gsl_multifit_fdfsolver_free (s);
}

void
test_ridge(void)
{
  const size_t n = 100;
  const size_t p = 10;
  const double xmin = -1.0;
  const double xmax = 1.0;
  const double dx = (xmax - xmin) / (n - 1.0);
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
  double *x = malloc(n * sizeof(double));
  double *y = malloc(n * sizeof(double));
  gsl_matrix *X = gsl_matrix_alloc(n, p);
  size_t i, j;

  /* construct artificial data */
  for (i = 0; i < n; ++i)
    {
      double ei = 0.2 * gsl_rng_uniform(r);

      x[i] = xmin + dx * i;
      y[i] = 1.0 / (1.0 + 25.0*x[i]*x[i]) + ei;
    }

  /* construct least squares matrix with polynomial model */
  for (i = 0; i < n; ++i)
    {
      double Xij = 1.0;

      for (j = 0; j < p; ++j)
        {
          gsl_matrix_set(X, i, j, Xij);
          Xij *= x[i];
        }
    }

  /* least squares fits */
  {
    gsl_multifit_linear_workspace *w = gsl_multifit_linear_alloc(n, p);
    gsl_vector_view yv = gsl_vector_view_array(y, n);
    gsl_vector *c0 = gsl_vector_alloc(p);
    gsl_vector *c1 = gsl_vector_alloc(p);
    gsl_vector *c2 = gsl_vector_alloc(p);
    gsl_vector *gamma = gsl_vector_calloc(p);
    gsl_matrix *cov = gsl_matrix_alloc(p, p);
    double chisq;

    /* test that ridge equals OLS solution for gamma = 0 */
    gsl_multifit_linear(X, &yv.vector, c0, cov, &chisq, w);
    gsl_multifit_linear_ridge(0.0, X, &yv.vector, c1, cov, &chisq, w);

    /* test c0 = c1 */
    for (j = 0; j < p; ++j)
      {
        double c0j = gsl_vector_get(c0, j);
        double c1j = gsl_vector_get(c1, j);

        gsl_test_rel(c1j, c0j, 1.0e-10, "test_ridge: gamma = 0, c0/c1");
      }

    for (i = 0; i < 7; ++i)
      {
        double g = pow(10.0, -(double) i);

        gsl_multifit_linear_ridge(g*g, X, &yv.vector, c1, cov,
                                  &chisq, w);

        gsl_vector_set_all(gamma, g);
        gsl_multifit_linear_ridge2(gamma, X, &yv.vector, c2, cov,
                                   &chisq, w);

        /* test c1 = c2 */
        for (j = 0; j < p; ++j)
          {
            double c1j = gsl_vector_get(c1, j);
            double c2j = gsl_vector_get(c2, j);

            gsl_test_rel(c2j, c1j, 1.0e-10, "test_ridge: gamma = %.1e", g);
          }
      }

    gsl_multifit_linear_free(w);
    gsl_vector_free(c0);
    gsl_vector_free(c1);
    gsl_vector_free(c2);
    gsl_vector_free(gamma);
    gsl_matrix_free(cov);
  }

  gsl_rng_free(r);
  free(x);
  free(y);
  gsl_matrix_free(X);
}
