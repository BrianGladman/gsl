/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <stdio.h>
#include <gsl_math.h>
#include <gsl_test.h>
#include <gsl_sum.h>

#define N 50

int
main (void)
{
  double t[N];
  double qnum[N];
  double qden[N];
  double dqnum[N*N];
  double dqden[N*N];
  double sum_accel;
  double dsum[N];
  double sum_plain;
  double prec;
  double x;
  double sd_actual, sd_est;
  int n;
  int i;

  for (i = 0; i < N*N; i++) 
    {
      dqnum[i] = (1+i) * 1e100 ;
      dqden[i] = (1+i) * 1e49 ;
    }

  {
    const double zeta_2 = M_PI * M_PI / 6.0;

    /* terms for zeta(2) */

    for (n = 0; n < N; n++)
      {
	double np1 = n + 1.0;
	t[n] = 1.0 / (np1 * np1);
      }

    gsl_sum_levin_u_accel (t, N, qnum, qden, &sum_accel, &sum_plain, &prec);
    gsl_test_rel (sum_accel, zeta_2, 1e-10, "Levin u-transform for zeta(2)");

    gsl_sum_levin_u_with_derivs_accel (t, N, qnum, qden, dqnum, dqden, dsum, 
				       &sum_accel, &sum_plain, &prec);
    gsl_test_rel (sum_accel, zeta_2, 1e-10, "transform with roundoff for zeta(2)");
    
    sd_est    = -log10(prec) ; 
    sd_actual = -log10(fabs(sum_accel - zeta_2)/zeta_2) ;

    gsl_test_abs (sd_est, sd_actual, 1.0, "error estimate (sd) for zeta(2)");
  }

  /* terms for exp(10.0) */
  x = 10.0;
  t[0] = 1.0;
  for (n = 1; n < N; n++)
    {
      t[n] = t[n - 1] * (x / n);
    }

  gsl_sum_levin_u_accel (t, N, qnum, qden, &sum_accel, &sum_plain, &prec);
  gsl_test_rel (sum_accel, exp (x), 1e-12, "Levin u-transform for exp(10.0)");

  gsl_sum_levin_u_with_derivs_accel (t, N, qnum, qden, dqnum, dqden, dsum, 
				     &sum_accel, &sum_plain, &prec);
  gsl_test_rel (sum_accel, exp (x), 1e-12, "Levin u-transform for exp(10.0)");


  /* terms for exp(-10.0) */
  x = -10.0;
  t[0] = 1.0;
  for (n = 1; n < N; n++)
    {
      t[n] = t[n - 1] * (x / n);
    }

  gsl_sum_levin_u_accel (t, N, qnum, qden, &sum_accel, &sum_plain, &prec);
  gsl_test_rel (sum_accel, exp (x), 1e-12, "Levin u-transform for exp(-10.0)");

  gsl_sum_levin_u_with_derivs_accel (t, N, qnum, qden, dqnum, dqden, dsum, 
				     &sum_accel, &sum_plain, &prec);
  gsl_test_rel (sum_accel, exp (x), 1e-12, "Levin u-transform for exp(-10.0)");

  sd_est    = -log10(prec) ; 
  sd_actual = -log10(fabs(sum_accel - exp(x))/exp(x)) ;
  
  gsl_test_abs (sd_est, sd_actual, 1.0, "error estimate (sd) for exp(-10.0)");


  /* terms for -log(1-x) */
  x = 0.5;
  t[0] = x;
  for (n = 1; n < N; n++)
    {
      t[n] = t[n - 1] * (x * n) / (n + 1.0);
    }

  gsl_sum_levin_u_accel (t, N, qnum, qden, &sum_accel, &sum_plain, &prec);
  gsl_test_rel (sum_accel, M_LN2, 1e-12, "Levin u-transform for -log(1/2)");

  /* terms for -log(1-x) */
  x = -1.0;
  t[0] = x;
  for (n = 1; n < N; n++)
    {
      t[n] = t[n - 1] * (x * n) / (n + 1.0);
    }

  gsl_sum_levin_u_accel (t, N, qnum, qden, &sum_accel, &sum_plain, &prec);
  gsl_test_rel (sum_accel, -M_LN2, 1e-12, "Levin u-transform for -log(2)");

  /* terms for an alternating asymptotic series */
  t[0] = 3.0 / (M_PI * M_PI);
  for (n = 1; n < N; n++)
    {
      t[n] = -t[n - 1] * (4.0 * (n + 1.0) - 1.0) / (M_PI * M_PI);
    }

  gsl_sum_levin_u_accel (t, N, qnum, qden, &sum_accel, &sum_plain, &prec);
  gsl_test_rel (sum_accel, 0.192594048773, 1e-10,
		"Levin u-transform for asymptotic alternating series");

  return gsl_test_summary ();
}






