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
  double square, difference;
  int i;
  
  sum = square = difference = 0;
  
  /* find the mean */
  the_mean = imean(array, size);

  /*  printf("ivariance: mean=%f\n", the_mean);*/
  
  /* find sum of the squares */
  for (i=0; i<size; i++){
    difference = array[i] - the_mean;
    square = difference * difference;
    sum = sum + square;
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

double iest_variance (int *array, int size)
{
  
  /* takes an array of integers and finds the variance */

  double sum, the_mean;
  double square, difference;
  int i;
  
  sum = square = difference = 0;
  
  /* find the mean */
  the_mean = imean(array, size);

  /*  printf("ivariance: mean=%f\n", the_mean);*/
  
  /* find sum of the squares */
  for (i=0; i<size; i++){
    difference = array[i] - the_mean;
    square = difference * difference;
    sum = sum + square;
  }

  /* find variance */
  return sum / (size-1);
}

double dest_variance (double *array, int size)
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
  return sum / (size -1);
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


double iest_sd (int *array, int size)
{

  /* finds the standard deviation of an array of integers */

  double variance, sd;

  variance = iest_variance(array, size);

  sd = sqrt(variance);

  return sd;
}

double dest_sd (double *array, int size)
{

  /* finds the standard deviation of an array of doubles */

  double variance, sd;

  variance = dest_variance(array, size);

  sd = sqrt(variance);

  return sd;
}

double iipvariance(int *array1, int *array2, int size1, int size2)
{
  /* Find the pooled variance of two integer arrays */

  double var1, var2, pooled_variance;

  pooled_variance = 0;

  /* find the variances */
  var1 = iest_variance(array1, size1);
  var2 = iest_variance(array2, size2);

  /* calculate the pooled variance */
  pooled_variance = (((size1 - 1)*var1)+((size2-1)*var2)) / (size1+size2-2);
  
  return pooled_variance;

}


double iittest (int *array1, int *array2, int size1, int size2)
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
  sd1 = iest_sd(array1, size1);
  sd2 = iest_sd(array2, size2);
  pv = iipvariance (array1, array2, size1, size2);
  
  /* calculate the t statistic */
  t = (mean1-mean2)/(sqrt(pv*((1.0/size1)+(1.0/size2))));
  
  return t;

}


/* double ddttest (int *array1, int *array2, int size1, int size2);*/
