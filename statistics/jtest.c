#include <stdio.h>

#include "gsl_bstats.h"

/* Test program for mean.c.  JimDavies 7.96 */

int 
main (void)
{
  
  int groupone[20];
  int grouptwo[20];
  double groupa[14];
  double groupb[14];
  double sd1, var1, pv, mean1, t;
  int maximum, minumum; 

  /* sample sets of integers */

  groupone[0] = 17 ;
  groupone[1] = 18 ;
  groupone[2] = 16 ;
  groupone[3] = 18 ;
  groupone[4] = 12 ;
  groupone[5] = 20 ;
  groupone[6] = 18 ;
  groupone[7] = 20 ;
  groupone[8] = 20 ;
  groupone[9] = 22 ;
  groupone[10] = 20;
  groupone[11] = 10;
  groupone[12] = 8 ;
  groupone[13] = 12;
  groupone[14] = 16;
  groupone[15] = 16;
  groupone[16] = 18;
  groupone[17] = 20;
  groupone[18] = 18;
  groupone[19] = 21;

  grouptwo[0]  = 19;
  grouptwo[1]  = 20;
  grouptwo[2]  = 22;
  grouptwo[3]  = 24;
  grouptwo[4]  = 10;
  grouptwo[5]  = 25;
  grouptwo[6]  = 20;
  grouptwo[7]  = 22;
  grouptwo[8]  = 21;
  grouptwo[9]  = 23;
  grouptwo[10] = 20;
  grouptwo[11] = 10;
  grouptwo[12] = 12;
  grouptwo[13] = 14;
  grouptwo[14] = 12;
  grouptwo[15] = 20;
  grouptwo[16] = 22;
  grouptwo[17] = 24;
  grouptwo[18] = 23;
  grouptwo[19] = 17;

  /* sample sets of doubles */

  groupa[0]  = .0421;
  groupa[1]  = .0941;
  groupa[2]  = .1064;
  groupa[3]  = .0242;
  groupa[4]  = .1331;
  groupa[5]  = .0773;
  groupa[6]  = .0243;
  groupa[7]  = .0815;
  groupa[8]  = .1186;
  groupa[9]  = .0356;
  groupa[10] = .0728;
  groupa[11] = .0999;
  groupa[12] = .0614;
  groupa[13] = .0479;
 
  groupb[0]  = .1081;
  groupb[1]  = .0986;
  groupb[2]  = .1566;
  groupb[3]  = .1961;
  groupb[4]  = .1125;
  groupb[5]  = .1942;
  groupb[6]  = .1079;
  groupb[7]  = .1021;
  groupb[8]  = .1583;
  groupb[9]  = .1673;
  groupb[10] = .1675;
  groupb[11] = .1856;
  groupb[12] = .1688;
  groupb[13] = .1512;
  
  /* test gsl_stats_imean */
  
  mean1 = gsl_stats_imean(groupone, 20);
  if (mean1 == 17){
    printf("success: gsl_stats_imean (%f)\n", mean1);
  }
  else {
    printf("ERROR: gsl_stats_imean (calculated %f)\n", mean1);
  }
  
  mean1 = 0;
  
  /* test gsl_stats_dmean */
  
  mean1 = gsl_stats_dmean(groupa, 14);
  if (mean1 > .072 && mean1 < .073){
    printf("success: gsl_stats_dmean (%f)\n", mean1);
  }
  else {
    printf("ERROR:gsl_stats_dmean (calculated %f)\n", mean1);
  }
  
  mean1 = 0;
  
  /* test gsl_stats_ivariance */
  
  var1 = gsl_stats_ivariance(groupone, 20);
  
  if (var1 > 13.6 && var1 < 13.8){
    printf("success: gsl_stats_ivariance (%f)\n", var1);
  }
  else{
    printf("ERROR: gsl_stats_ivariance (calculated %f)\n", var1);
  }
  
  var1 = 0;
  
  /* test gsl_stats_dvariance */
  
  var1 = gsl_stats_dvariance(groupa, 14);
  
  if (var1 > .001135 && var1 < .001139){
    printf("success: gsl_stats_dvariance (%f)\n", var1);
  }
  else{
    printf("ERROR: gsl_stats_dvariance (calculated %f)\n", var1);
  }
  
  var1 = 0;
  
  /* test gsl_stats_iest_variance */
  
  var1 = gsl_stats_iest_variance(groupone, 20);
  
  if (var1 > 14.42 && var1 < 14.425){
    printf("success: gsl_stats_iest_variance (%f)\n", var1);
  }
  else{
    printf("ERROR: gsl_stats_iest_variance (calculated %f)\n", var1);
  }
  
  var1 = 0;
  
  /* test gsl_stats_dest_variance */
  
  var1 = gsl_stats_dest_variance(groupb, 14);
  
  if (var1 > .001 && var1 < .0013){
    printf("success: gsl_stats_dest_variance (%f)\n", var1);
  }
  else{
    printf("ERROR: gsl_stats_dest_variance (calculated %f)\n", var1);
  }
  
  var1 = 0;
  
  /* test gsl_stats_isd */
  
  sd1 = gsl_stats_isd(groupone, 20);
  if (sd1 > 3.6 && sd1 < 3.78){
    printf("success: gsl_stats_isd (%f)\n", sd1);
  }
  else{
    printf("ERROR: gsl_stats_isd (calculated %f)\n", sd1);
  }
  
  sd1 = 0;

  /* test gsl_stats_dsd */
  
  sd1 = gsl_stats_dsd(groupa, 14);
  if (sd1 > 0.0336 && sd1 < 0.0338){
    printf("success: gsl_stats_dsd (%f)\n", sd1);
  }
  else{
    printf("ERROR: gsl_stats_dsd (calculated %f)\n", sd1);
  }
  
  sd1 = 0;
  
  /* test gsl_stats_iest_sd */
  
  sd1 = gsl_stats_iest_sd(groupone, 20);
  if (sd1 > 3.79 && sd1 < 3.8){
    printf("success: gsl_stats_iest_sd (%f)\n", sd1);
  }
  else{
    printf("ERROR: gsl_stats_iest_sd (calculated %f)\n", sd1);
  }
  
  sd1 = 0;
  
  /* test gsl_stats_dest_sd */
  
  sd1 = gsl_stats_dest_sd(groupa, 14);
  if (sd1 > .034  && sd1 < .036){
    printf("success: gsl_stats_dest_sd (%f)\n", sd1);
  }
  else{
    printf("ERROR: gsl_stats_dest_sd (calculated %f)\n", sd1);
  }
  
  sd1 = 0;
  
  /* test gsl_stats_ipvariance */
  
  pv = gsl_stats_ipvariance(groupone, grouptwo, 20, 20);
  if (pv > 18.84 && pv < 18.85 ){
    printf("success: gsl_stats_ipvariance, (%f)\n", pv);
  }
  else{
    printf("ERROR: gsl_stats_ipvariance (calculated %f)\n", pv);
  }
  
  pv = 0;
  
  /* test gsl_stats_dpvariance */
  
  pv = gsl_stats_dpvariance(groupa, groupb, 14, 14);
  if (pv > 0.00122 && pv < 0.00124 ){
    printf("success: gsl_stats_dpvariance, (%f)\n", pv);
  }
  else
    printf("ERROR: gsl_stats_dpvariance (calculated %f)\n", pv);
  
  pv = 0;
  
  /* test gsl_stats_ittest */

  t = gsl_stats_ittest(groupone, grouptwo, 20, 20);
  if (t > -1.47 && t < -1.45){
    printf("success: gsl_stats_ittest (%f)\n", t);
  }
  else{
    printf("ERROR: gsl_stats_ittest (calculated %f)\n", t);
  }
  
  t = 0;
  
  /* test gsl_stats_dttest */
  
  t = gsl_stats_dttest(groupa, groupb, 14, 14);
  if (t > -5.68 && t < -5.66){
    printf("success: gsl_stats_dttest (%f)\n", t);
  }
  else{
    printf("ERROR: gsl_stats_dttest (calculated %f)\n", t);
  }
  
  t = 0;

  /* test gsl_stats_imax */
  
  maximum = gsl_stats_imax(groupone, 20);
  if (maximum == 22){
    printf("success: gsl_stats_imax (%d)\n", maximum);
  }
  else{
    printf("ERROR: gsl_stats_imax (calculated %d)\n", maximum);
  }
  
  maximum = 0;

  /* test gsl_stats_dmax */
  
  t = gsl_stats_dmax(groupa, 14);
  if (t > 0.132 && t < 0.134){
    printf("success: gsl_stats_dmax (%f)\n", t);
  }
  else{
    printf("ERROR: gsl_stats_dmax (calculated %f)\n", t);
  }
  
  t = 0;

  /* test gsl_stats_imin */
  
  minumum = gsl_stats_imin(groupone, 20);
  if (minumum == 8){
    printf("success: gsl_stats_imin (%d)\n", minumum);
  }
  else{
    printf("ERROR: gsl_stats_imin (calculated %d)\n", minumum);
  }
  
  minumum = 0;

  /* test gsl_stats_dmin */
  
  t = gsl_stats_dmin(groupa, 14);
  if (t > 0.023 && t < 0.025){
    printf("success: gsl_stats_dmin (%f)\n", t);
  }
  else{
    printf("ERROR: gsl_stats_dmin (calculated %f)\n", t);
  }
  
  t = 0;



  return 0;
}  






