#include <stdio.h>
#include <math.h>

#include "variance.h"
#include "mean.h"


double gsl_stats_ivariance (int *array, int size)
{
  
  /* takes an array of integers and finds the variance */

  double sum, the_mean;
  double square, difference;
  int i;
  
  sum = square = difference = 0;
  
  /* find the mean */
  the_mean = gsl_stats_imean(array, size);

  /* find sum of the squares */
  for (i=0; i<size; i++){
    difference = array[i] - the_mean;
    square = difference * difference;
    sum = sum + square;
  }

  /* find variance */
  return sum / size;
}

double gsl_stats_dvariance (double *array, int size)
{
  
  /* takes an array of doubles and finds the variance */

  double sum, the_mean;
  int i;

  sum=0;
  
  /* find the mean */
  the_mean = gsl_stats_dmean(array, size);

  /* find sum of the squares */
  for (i=0; i<size; i++){
    sum += ((array[i] - the_mean) * (array[i] - the_mean));
  }

  /* find variance */
  return (sum / size);
}

double gsl_stats_iest_variance (int *array, int size)
{
  /* takes an array of integers and finds the variance */

  double sum, the_mean;
  double square, difference;
  int i;
  
  sum = square = difference = 0;
  
  /* find the mean */
  the_mean = gsl_stats_imean(array, size);

  /* find sum of the squares */
  for (i=0; i<size; i++){
    difference = array[i] - the_mean;
    square = difference * difference;
    sum = sum + square;
  }

  /* find variance */
  return sum / (size-1);
}

double gsl_stats_dest_variance (double *array, int size)
{
  /* takes an array of doubles and finds the variance */

  double sum, the_mean;
  int i;
  
  sum=0;
  
  /* find the mean */
  the_mean = gsl_stats_dmean(array, size);
  
  /* find sum of the squares */
  for (i=0; i<size; i++){
    sum += ((array[i] - the_mean) * (array[i] - the_mean));
  }

  /* find variance */
  return sum / (size -1);
}

double gsl_stats_isd (int *array, int size)
{
  /* finds the standard deviation of an array of integers */

  double variance, sd;

  variance = gsl_stats_ivariance(array, size);

  sd = sqrt(variance);

  return sd;
}

double gsl_stats_dsd (double *array, int size)
{
  /* finds the standard deviation of an array of doubles */

  double variance, sd;

  variance = gsl_stats_dvariance(array, size);

  sd = sqrt(variance);

  return sd;
}

double gsl_stats_iest_sd (int *array, int size)
{
  /* finds the standard deviation of an array of integers */

  double variance, sd;

  variance = gsl_stats_iest_variance(array, size);

  sd = sqrt(variance);

  return sd;
}

double gsl_stats_dest_sd (double *array, int size)
{
  /* finds the standard deviation of an array of doubles */

  double variance, sd;

  variance = gsl_stats_dest_variance(array, size);

  sd = sqrt(variance);

  return sd;
}

double gsl_stats_iipvariance(int *array1, int *array2, int size1, int size2)
{
  /* Find the pooled variance of two integer arrays */

  double var1, var2, pooled_variance;

  pooled_variance = 0;

  /* find the variances */
  var1 = gsl_stats_iest_variance(array1, size1);
  var2 = gsl_stats_iest_variance(array2, size2);

  /* calculate the pooled variance */
  pooled_variance = (((size1 - 1)*var1)+((size2-1)*var2)) / (size1+size2-2);
  
  return pooled_variance;

}

double gsl_stats_ddpvariance(double *array1, double *array2, int size1, int size2)
{
  /* Find the pooled variance of two double arrays */

  double var1, var2, pooled_variance;

  pooled_variance = 0;

  /* find the variances */
  var1 = gsl_stats_dest_variance(array1, size1);
  var2 = gsl_stats_dest_variance(array2, size2);

  /* calculate the pooled variance */
  pooled_variance = (((size1 - 1)*var1)+((size2-1)*var2)) / (size1+size2-2);

  return pooled_variance;

}


int gsl_stats_imax (int *array, int size)
{
  /* finds the highest member of an integer array */
  int max, i;
  
  max = array[0];
  
  for (i=0; i < size; i++){
    if (array[i] > max) 
	      max = array[i];
  }
  return max;
}

double gsl_stats_dmax (double *array, int size)
{
  /* finds the highest member of an integer array */
  double max;
  int i;

  max = array[0];
  
  for (i=0; i < size; i++){
    if (array[i] > max)
	      max = array[i];
  }
  return max;
}

int gsl_stats_imin (int *array, int size)
{
    /* finds the highest member of an integer array */
  
  int min, i;

  min = array[0];
  
  for (i=0; i < size; i++){
    if (array[i] < min)
	      min = array[i];
  }
  return min;
}

double gsl_stats_dmin (double *array, int size)
{
  /* finds the highest member of an integer array */
  
  double min;
  int i;
  
  min = array[0];
  
  for (i=0; i < size; i++){
    if (array[i] <min)
      min = array[i];
  }

  return min;
}

