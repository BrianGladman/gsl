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

double fmean (float *array, int size)
{
  
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
  /*printf("MEAN.C: The standard deviation is: %f.\n", sd);*/

  return sd;
}

double dsd (double *array, int size)
{

  /* finds the standard deviation of an array of doubles */

  double variance, sd;

  variance = dvariance(array, size);

  sd = sqrt(variance);
  printf("MEAN.C: The standard deviation is: %f.\n", sd);

  return sd;
}
