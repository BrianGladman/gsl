#include <stdio.h>

#include "mean.h"

/* Test program for mean.c.  JimDavies 7.96 */

main()
{
  
  int grades[5];
  double measure[6];
  double average, var, sd_grades;
  double mean_measure, var_measure, sd_measure;

  /* integers test */

  grades[0] = 99;
  grades[1] = 76;
  grades[2] = 95;
  grades[3] = 85;
  grades[4] = 457;
   
  average = imean(grades, 5);
  /*   printf ("Your average is %f.\n", average); */
  var = ivariance(grades, 5);
  /*  printf ("The variance of the grades is %f.\n\n", var);*/
  sd_grades = isd(grades, 5);
  /*  printf("JTEST: the SD is: %f.\n", sd_grades);*/

  /* doubles test */

  measure[0] = 24.655;
  measure[1] = 20.1;
  measure[2] = 30.001;
  measure[3] = 25;
  measure[4] = 24.99;
  measure[5] = 15.67;

  mean_measure = dmean(measure, 6);
  var_measure = dvariance(measure, 6);
  sd_measure = dsd(measure, 6);
  printf("jtest: the sd is: %f.\n", sd_measure);
  
  return 0;
}
