#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_interp.h>

typedef double (*eval_func)(const gsl_interp*,const double *x,const double *y,
                            double,gsl_interp_accel*);

static const double EPSILON = 0.00000001;
static const double EPS = 0.0001;
static const int PTS = 801;

static const double x_values[] = { 0, 1, 2, 3, 4, 5 };
static const double y_values[] = { 5, 0, 5, 2, 0, 5 };
static double x_range;

static const int N = sizeof(x_values)/sizeof(double);

static double* interpolate(const gsl_interp_type *type, const size_t order);
static void    check_discont(gsl_interp *interp, const eval_func func,
                             const char *dsym, const double x,
                             gsl_interp_accel *accel);

static inline double
scale_x(const int i)
{
  return x_values[0] + (double)i/(PTS-1)*x_range;
}

int
main(int argc, char *argv[])
{
  int i;
  double *cspline, *cspline_d, *cspline_dd;
  double *akima, *akima_d, *akima_dd;

  x_range = x_values[N-1] - x_values[0];

  akima = interpolate(gsl_interp_akima_periodic, 0);
  akima_d = interpolate(gsl_interp_akima_periodic, 1);
  akima_dd = interpolate(gsl_interp_akima_periodic, 2);
  cspline = interpolate(gsl_interp_cspline_periodic, 0);
  cspline_d = interpolate(gsl_interp_cspline_periodic, 1);
  cspline_dd = interpolate(gsl_interp_cspline_periodic, 2);

  puts("# x   akima  cspline  akima'  cspline'  akima''  cspline''");
  for (i = -PTS/2; i < PTS/2; i++)
    printf("%0.5f %0.5f %0.5f %0.5f %0.5f %0.5f %0.5f\n",
           scale_x(i),
           akima[(i + PTS)%PTS],
           cspline[(i + PTS)%PTS],
           akima_d[(i + PTS)%PTS],
           cspline_d[(i + PTS)%PTS],
           akima_dd[(i + PTS)%PTS],
           cspline_dd[(i + PTS)%PTS]);

  free(akima);
  free(akima_d);
  free(akima_dd);
  free(cspline);
  free(cspline_d);
  free(cspline_dd);
  return 0;
}

static double*
interpolate(const gsl_interp_type *type, const size_t order)
{
  static const eval_func DERIV_ORDERS[] = {
    gsl_interp_eval,
    gsl_interp_eval_deriv,
    gsl_interp_eval_deriv2
  };
  static const char *DERIV_SYMS[] = {
    "",
    "'",
    "''"
  };

  gsl_interp *interp = gsl_interp_alloc(type, N);
  gsl_interp_accel *accel = gsl_interp_accel_alloc();
  double *plotpts = (double*)malloc(PTS*sizeof(double));
  int i;

  assert(order < sizeof(DERIV_ORDERS));
  gsl_interp_init(interp, x_values, y_values, N);
  for (i = 0; i < PTS; i++)
    plotpts[i] = DERIV_ORDERS[order](interp, x_values, y_values,
                                     scale_x(i), accel);

  for (i = 0; i < N-1; i++)
    check_discont(interp, DERIV_ORDERS[order], DERIV_SYMS[order],
                  x_values[i], accel);

  gsl_interp_free(interp);
  gsl_interp_accel_free(accel);

  return plotpts;
}

static void
check_discont(gsl_interp *interp, const eval_func func, const char *dsym,
              const double x, gsl_interp_accel *accel)
{
  double xm, xp;
  double ym, yp;

  xm = fmod(x - EPSILON + x_range, x_range);
  xp = fmod(x + EPSILON, x_range);
  ym = func(interp, x_values, y_values, xm, accel);
  yp = func(interp, x_values, y_values, xp, accel);
  if (fabs(yp - ym) > EPS) {
    fprintf(stderr, "discontinuity in %s%s at %f:\n"
                    "\t left value: %f\n"
                    "\tright value: %f\n"
                    "\t difference: %f\n\n",
            gsl_interp_name(interp), dsym, x, ym, yp, yp-ym);
  }
}
