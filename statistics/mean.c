#include <stdio.h>
#include <math.h>

#include "mean.h"

double imean (int *array, int size)
{

  /* Takes an array of integers and finds a double mean */
  
  double sum, mean;
  int i;
  
  sum = 0;
  
  /*  Find the sum    */
  for (i=0; i<size; i++){
    sum += array[i];
  }
  
  /* find the mean */
  mean = (sum / size);
  
  return mean;
}

double dmean (double *array, int size)
{

  /* Takes an array of doubles and finds the double mean */
  
  double sum, mean;
  int i;
  
  sum = 0;
  
  /*  Find the sum    */
  for (i=0; i<size; i++){
    sum += array[i];
  }
  
  /* find the mean */
  mean = (sum / size);

  return mean;
}

double ivariance (int *array, int size)
{
  
  /* takes an array of integers and finds the variance */

  double sum, the_mean;
  int i;
  
  sum=0;
  
  /* find the mean */
  the_mean = imean(array, size);
  
  /* find sum of the squares */
  for (i=0; i<size; i++){
    sum += ((array[i] - the_mean) * (array[i] - the_mean));
  }

  /* find variance */
  return sum / size;
}

double dvariance (double *array, int size)
{
  
  /* takes an array of doubles and finds the variance */

  double sum, the_mean;
  int i;
  
  sum=0;
  
  /* find the mean */
  the_mean = dmean(array, size);
  
  /* find sum of the squares */
  for (i=0; i<size; i++){
    sum += ((array[i] - the_mean) * (array[i] - the_mean));
  }

    /* find variance */
  return (sum / size);
}

double isd (int *array, int size)
{

  /* finds the standard deviation of an array of integers */

  double variance, sd;

  variance = ivariance(array, size);

  sd = sqrt(variance);

  return sd;
}

double dsd (double *array, int size)
{

  /* finds the standard deviation of an array of doubles */

  double variance, sd;

  variance = dvariance(array, size);

  sd = sqrt(variance);

  return sd;
}

double iipvariance(int *array1, int *array2, int size1, int size2)
{
  /* Find the pooled variance of two integer arrays */

  double var1, var2, pooled_variance;
  
  /* find the variances */
  var1 = ivariance(array1, size1);
  var2 = ivariance(array2, size2);

  /* calculate the pooled variance */
  pooled_variance = (((size1 - 1)*var1)+((size2-1)*var2)) / (sqrt(size1+size2-2));
  
  return pooled_variance;

}


double iittest (int *array1, int *array2, int size1, int size2, double alpha, int tails)
{

  /* runs a t-test between two arrays of integers representing 
      independent samples. Tests to see if the difference between
      means of the samples is different from zero */

  double mean1, mean2;  /* means of the two samples */
  double sd1, sd2;  /* standard deviations */
  double pv;        /* pooled variance */
  double t;         /* the t statistic */
  
  /* find means and standard deviations for the two samples */
  mean1 = imean(array1, size1);
  mean2 = imean(array2, size2);
  sd1 = isd(array1, size1);
  sd2 = isd(array2, size2);
  pv = iipvariance (array1, array2, size1, size2);
  
  /* calculate the t statistic */
  t = (mean1-mean2)/(sqrt(pv*((1/size1)+(1/size2))));

}
