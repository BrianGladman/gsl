/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <stdio.h>
#include <gsl_math.h>
#include <gsl_test.h>
#include <gsl_sum.h>

#define N 50

void check_trunc (double * t, double expected, const char * desc);
void check_full (double * t, double expected, const char * desc);

int
main (void)
{
  {
    double t[N];
    int n;

    const double zeta_2 = M_PI * M_PI / 6.0;

    /* terms for zeta(2) */

    for (n = 0; n < N; n++)
      {
	double np1 = n + 1.0;
	t[n] = 1.0 / (np1 * np1);
      }

    check_trunc (t, zeta_2, "zeta(2)");
    check_full (t, zeta_2, "zeta(2)");
  }

  {
    double t[N];
    double x, y;
    int n;

    /* terms for exp(10.0) */
    x = 10.0;
    y = exp(x);

    t[0] = 1.0;
    for (n = 1; n < N; n++)
      {
	t[n] = t[n - 1] * (x / n);
      }

    check_trunc (t, y, "exp(10)");
    check_full (t, y, "exp(10)");
  }

  {
    double t[N];
    double x, y;
    int n;

    /* terms for exp(-10.0) */
    x = -10.0;
    y = exp(x);

    t[0] = 1.0;
    for (n = 1; n < N; n++)
      {
	t[n] = t[n - 1] * (x / n);
      }

    check_trunc (t, y, "exp(-10)");
    check_full (t, y, "exp(-10)");
  }

  {
    double t[N];
    double x, y;
    int n;

    /* terms for -log(1-x) */
    x = 0.5;
    y = -log(1-x);
    t[0] = x;
    for (n = 1; n < N; n++)
      {
	t[n] = t[n - 1] * (x * n) / (n + 1.0);
      }

    check_trunc (t, y, "-log(1/2)");
    check_full (t, y, "-log(1/2)");
  }

  {
    double t[N];
    double x, y;
    int n;

    /* terms for -log(1-x) */
    x = -1.0;
    y = -log(1-x);
    t[0] = x;
    for (n = 1; n < N; n++)
      {
	t[n] = t[n - 1] * (x * n) / (n + 1.0);
      }

    check_trunc (t, y, "-log(2)");
    check_full (t, y, "-log(2)");
  }

  {
    double t[N];
    int n;

    double result = 0.192594048773;

    /* terms for an alternating asymptotic series */

    t[0] = 3.0 / (M_PI * M_PI);

    for (n = 1; n < N; n++)
      {
	t[n] = -t[n - 1] * (4.0 * (n + 1.0) - 1.0) / (M_PI * M_PI);
      }

    check_trunc (t, result, "asymptotic series");
    check_full (t, result, "asymptotic series");
  }

  return gsl_test_summary ();
}

void
check_trunc (double * t, double expected, const char * desc)
{
  double qnum[N], qden[N];
  double sum_accel, sum_plain, prec, sd_actual, sd_est;
  size_t n_used;
  
  gsl_sum_levin_u_trunc_accel (t, N, qnum, qden,
                               &sum_accel, &n_used, &sum_plain, &prec);
  gsl_test_rel (sum_accel, expected, 1e-10, "trunc result, %s", desc);
  
  sd_est = -log10 (prec);
  sd_actual = -log10 (DBL_EPSILON + fabs (sum_accel - expected) / expected);
  gsl_test (sd_est > sd_actual, "trunc signficant digits, %s", desc);
}

void
check_full (double * t, double expected, const char * desc)
{
  double qnum[N], qden[N], dqnum[N * N], dqden[N * N], dsum[N];
  double sum_accel, sum_plain, prec, sd_actual, sd_est;
  size_t n_used;
  
  gsl_sum_levin_u_accel (t, N, qnum, qden, dqnum, dqden, dsum,
                         &sum_accel, &n_used, &sum_plain, &prec);
  gsl_test_rel (sum_accel, expected, 1e-10, "full result, %s", desc);
  
  sd_est = -log10 (prec);
  sd_actual = -log10 (DBL_EPSILON + fabs (sum_accel - expected) / expected);
  gsl_test (sd_est > sd_actual, "full significant digits, %s", desc);
}
