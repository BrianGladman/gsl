#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sum.h>

#define N 20

int 
main (void)
{
  double t[N], qnum[N], qden[N], dsum[N], dqnum[N * N], dqden[N * N];
  double sum_accel, sum_plain, prec;
  double sum = 0;
  size_t n_used ;
  int n;
  
  const double zeta_2 = M_PI * M_PI / 6.0;
  
  /* terms for zeta(2) = \sum_{n=0}^{\infty} 1/n^2 */

  for (n = 0; n < N; n++)
    {
      double np1 = n + 1.0;
      t[n] = 1.0 / (np1 * np1);
      sum += t[n] ;
    }
  
  gsl_sum_levin_u_accel (t, N, qnum, qden, dqnum, dqden, dsum,
                         &sum_accel, &n_used, &sum_plain, &prec);


  printf("term-by-term sum = % .16f using %d terms\n", sum, N) ;
  printf("term-by-term sum = % .16f using %d terms\n", sum_plain, n_used) ;
  printf("exact value      = % .16f\n", zeta_2) ;
  printf("accelerated sum  = % .16f using %d terms\n", sum_accel, n_used) ;
  printf("estimated error  = % .16f\n", prec * sum_accel) ;
  printf("actual error     = %+.16f\n", sum_accel - zeta_2) ;

        
  return 0;
}
