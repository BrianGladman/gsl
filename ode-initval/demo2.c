#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

gsl_odeiv_control *
gsl_odeiv_control_scaled_new(double eps_abs, double eps_rel,
                             double a_y, double a_dydt,
                             const double scale_y[],
                             size_t dim);


int
func (double t, const double y[], double f[], void *params)
{
  double mu = *(double *)params;
  f[0] = y[1];
  f[1] = -y[0] - mu*y[1]*(y[0]*y[0] - 1);
  return GSL_SUCCESS;
}

int
jac (double t, const double y[], double *dfdy, double dfdt[],
            void *params)
{
  double mu = *(double *)params;
  gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, 2, 2);
  gsl_matrix_set (&dfdy_mat.matrix, 0, 0, 0.0);
  gsl_matrix_set (&dfdy_mat.matrix, 0, 1, 1.0);
  gsl_matrix_set (&dfdy_mat.matrix, 1, 0, -2.0*mu*y[0]*y[1] - 1.0);
  gsl_matrix_set (&dfdy_mat.matrix, 1, 1, -mu*(y[0]*y[0] - 1.0));
  dfdt[0] = 0.0;
  dfdt[1] = 0.0;
  return GSL_SUCCESS;
}

int
main (void)
{
  const gsl_odeiv_step_type * T = gsl_odeiv_step_rk2;

  gsl_odeiv_step * s = gsl_odeiv_step_alloc (T, 2);
  
  double scale_y[2] = { 1e-12, 2.0 };

  gsl_odeiv_control * c = gsl_odeiv_control_scaled_new (1e-6, 1e-6,
                                                        1, 0,
                                                        scale_y,
                                                        2);
  gsl_odeiv_evolve * e = gsl_odeiv_evolve_alloc (2);
  double mu = 10;
  gsl_odeiv_system sys = {func, jac, 2, &mu};

  int i;
  double t = 0.0, t1 = 100.0;
  double h = 1e-6;
  double y[2] = { 1.0, 0.0 } ;

  gsl_ieee_env_setup();

#if 0
  while (t < t1)
    {
      int status = gsl_odeiv_evolve_apply(e, c, s, &sys, &t, t1, &h, y);

      if (status != GSL_SUCCESS)
        {
          break;
        }

      printf("%.5e %.5e %.5e\n", t, y[0], y[1]);
    }
#else
  for (i = 1; i <= 100 ; i++)
    {
      double ti = i * t1 / 100.0;

        while (t < ti)
          {
            int status = gsl_odeiv_evolve_apply(e, c, s, &sys, &t, ti, &h, y);
            
            if (status != GSL_SUCCESS)
              {
                break;
              }
          }
      printf("%.5e %.5e %.5e\n", t, y[0], y[1]);
    }
#endif

  gsl_odeiv_evolve_free(e);
  gsl_odeiv_control_free(c);
  gsl_odeiv_step_free(s);
}
