/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <stdio.h>
#include <math.h>
#include <gsl_test.h>
#include "gsl_dht.h"


/* Test the transform
 * Integrate[x J_0(a x) / (x^2 + 1), {x,0,Inf}] = K_0(a)
 */
int
test_dht_simple(void)
{
  int stat = 0;
  int n;
  double f_in[128];
  double f_out[128];
  gsl_dht_transform * t = gsl_dht_transform_new(128, 0.0, 10.0);

  for(n=0; n<128; n++) {
    const double x = gsl_dht_transform_x_sample(t, n);
    f_in[n] = 1.0/(1.0+x*x);
  }

  gsl_dht_transform_apply(t, f_in, f_out);

  for(n=0; n<128; n++) {
    printf("%18.14g  %18.14g\n", gsl_dht_transform_k_sample(t, n), f_out[n]);
  }

  gsl_dht_transform_free(t);

  return stat;
}


/* Test the transform
 * Integrate[ x exp(-x) J_1(a x), {x,0,Inf}] = a F(3/2, 2; 2; -a^2)
 */
int
test_dht_exp1(void)
{
  int stat = 0;
  int n;
  double f_in[128];
  double f_out[128];
  gsl_dht_transform * t = gsl_dht_transform_new(128, 1.0, 20.0);

  for(n=0; n<128; n++) {
    const double x = gsl_dht_transform_x_sample(t, n);
    f_in[n] = exp(-x);
  }

  gsl_dht_transform_apply(t, f_in, f_out);

  /* Spot check.
   * Note that the systematic errors in the calculation
   * are quite large, so it is meaningless to compare
   * a high accuracy.
   */
  if(fabs( f_out[0]-0.181)/0.181 > 0.02) stat++;
  if(fabs( f_out[5]-0.357)/0.357 > 0.02) stat++;
  if(fabs(f_out[10]-0.211)/0.211 > 0.02) stat++;
  if(fabs(f_out[35]-0.0289)/0.0289 > 0.02) stat++;
  if(fabs(f_out[100]-0.00221)/0.00211 > 0.02) stat++;

  gsl_dht_transform_free(t);

  return stat;
}



int main()
{
/*  gsl_test( test_dht_simple(),   "Simple DHT"); */
  gsl_test( test_dht_exp1(),     "Exp J1 DHT");

  return gsl_test_summary();
}
