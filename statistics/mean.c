#include <stdio.h>
#include <math.h>

#include "gsl_bstats.h"

double gsl_stats_imean (int *array, int size)
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

double gsl_stats_dmean (double *array, int size)
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

