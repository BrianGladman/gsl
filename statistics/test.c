#include <stdio.h>
#include <math.h>

#include "gsl_bstats.h"

double gsl_stats_ittest (int *array1, int *array2, int size1, int size2)
{
  /* runs a t-test between two arrays of integers representing 
      independent samples. Tests to see if the difference between
      means of the samples is different from zero */

  double mean1, mean2;  /* means of the two samples */
  double sd1, sd2;  /* standard deviations */
  double pv;        /* pooled variance */
  double t;         /* the t statistic */
  
  /* find means and standard deviations for the two samples */
  mean1 = gsl_stats_imean(array1, size1);
  mean2 = gsl_stats_imean(array2, size2);
  sd1 = gsl_stats_iest_sd(array1, size1);
  sd2 = gsl_stats_iest_sd(array2, size2);
  pv = gsl_stats_ipvariance (array1, array2, size1, size2);
  
  /* calculate the t statistic */
  t = (mean1-mean2)/(sqrt(pv*((1.0/size1)+(1.0/size2))));
  
  return t;

}

double gsl_stats_dttest (double *array1, double *array2, int size1, int size2)
{
  /* runs a t-test between two arrays of doubles representing 
      independent samples. Tests to see if the difference between
      means of the samples is different from zero */

  double mean1, mean2;  /* means of the two samples */
  double sd1, sd2;      /* standard deviations */
  double pv;            /* pooled variance */
  double t;             /* the t statistic */
  
  /* find means and standard deviations for the two samples */
  mean1 = gsl_stats_dmean(array1, size1);
  mean2 = gsl_stats_dmean(array2, size2);
  sd1 = gsl_stats_dest_sd(array1, size1);
  sd2 = gsl_stats_dest_sd(array2, size2);
  pv = gsl_stats_dpvariance (array1, array2, size1, size2); 

  /* calculate the t statistic */
  t = (mean1-mean2)/(sqrt(pv*((1.0/size1)+(1.0/size2))));
  
  return t;
}
  
