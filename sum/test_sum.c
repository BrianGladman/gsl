/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <stdio.h>
#include <gsl_math.h>
#include <gsl_test.h>
#include <gsl_sum.h>


double frac_diff(double x1, double x2)
{
  if(x1 <= DBL_MAX && x2 <= DBL_MAX)
    return fabs((x1-x2)/(x1+x2));
  else
    return 1.0;
}

#define N 35

int check_levin_u(void)
{
  const double zeta_2 = M_PI*M_PI/6.0;
  double array[N];
  double qnum[N];
  double qden[N];
  double sum_accel;
  double sum_plain;
  double prec;
  double x;
  int n;

  int status = 0;
  int s;

  /* terms for zeta(2) */
  for(n=0; n<N; n++) {
    double np1 = n + 1.0;
    array[n] = 1.0/(np1*np1);
  }

  s = 0;
  gsl_sum_levin_u_impl(array, N, qnum, qden, &sum_accel, &sum_plain, &prec);
  s += ( frac_diff( sum_accel, zeta_2 ) > 1.0e-10 );
  gsl_test(s, "  Levin u-Transform for zeta(2)");
  status += s;

  /* terms for exp(10.0) */
  x = 10.0;
  array[0] = 1.0;
  for(n=1; n<N; n++) {
    array[n] = array[n-1] * (x / n);
  }

  s = 0;
  gsl_sum_levin_u_impl(array, N, qnum, qden, &sum_accel, &sum_plain, &prec);
  s += ( frac_diff( sum_accel, exp(x)) > 1.0e-12 );
  gsl_test(s, "  Levin u-Transform for exp(10.0)");
  status += s;

  /* terms for exp(-10.0) */
  x = -10.0;
  array[0] = 1.0;
  for(n=1; n<N; n++) {
    array[n] = array[n-1] * (x / n);
  }

  s = 0;
  gsl_sum_levin_u_impl(array, N, qnum, qden, &sum_accel, &sum_plain, &prec);
  s += ( frac_diff( sum_accel, exp(x)) > 1.0e-12 );
  printf("%.18g %g\n%.18g\n%.18g\n", sum_accel, prec, sum_plain, exp(x)) ;
  gsl_test(s, "  Levin u-Transform for exp(-10.0)");
  status += s;

  /* terms for -log(1-x) */
  x = 0.5;
  array[0] = x;
  for(n=1; n<N; n++) {
    array[n]  = array[n-1] * (x * n) / (n+1.0) ;
  }

  s = 0;
  gsl_sum_levin_u_impl(array, N, qnum, qden, &sum_accel, &sum_plain, &prec);
  s += ( frac_diff( sum_accel, M_LN2) > 1.0e-12 );
  gsl_test(s, "  Levin u-Transform for -log(1/2)");
  status += s;


  /* terms for -log(1-x) */
  x = -1.0;
  array[0] = x;
  for(n=1; n<N; n++) {
    array[n]  = array[n-1] * (x * n) / (n+1.0) ;
  }

  s = 0;
  gsl_sum_levin_u_impl(array, N, qnum, qden, &sum_accel, &sum_plain, &prec);
  s += ( frac_diff( sum_accel, -M_LN2) > 1.0e-12 );
  gsl_test(s, "  Levin u-Transform for -log(2)");
  status += s;


  /* terms for an alternating asymptotic series */
  array[0] = 3.0/(M_PI*M_PI);
  for(n=1; n<N; n++) {
    array[n] = -array[n-1] * (4.0*(n+1.0) - 1.0)/(M_PI*M_PI);
  }

  s = 0;
  gsl_sum_levin_u_impl(array, N, qnum, qden, &sum_accel, &sum_plain, &prec);
  printf("%.18g %g\n%.18g\n%.18g\n", sum_accel, prec, sum_plain,0.192594048773) ;
  s += ( frac_diff( sum_accel, 0.192594048773) > 1.0e-10 );
  gsl_test(s, "  Levin u-Transform for asymptotic alternating series");
  status += s;


  return status;
}


int main(int argc, char * argv[])
{
  argc = 0 ; /* prevent warnings about unused parameters */
  argv = 0 ;

  gsl_test(check_levin_u(),      "Levin u-Transform");

  return gsl_test_summary();
}
