#include <stdio.h>

#include "mean.h"

/* Test program for mean.c.  JimDavies 7.96 */

main()
{
  
  int grades[5];
  int groupone[20];
  int grouptwo[20];
  double measure[6];
  double average, var, sd_grades;
  double mean_measure, var_measure, sd_measure;
  double sd1, var1, pv, mean1;
  
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
  /*   printf("jtest: the sd is: %f.\n", sd_measure);*/

  /* integers ttest */

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
 
  /* test imean */
 
  mean1 = imean(groupone, 20);
  if (mean1 = 17){
    printf("success: imean\n");
  }
  else {
    printf("ERROR: imean (calculated %f)\n", mean1);
  }
  
  mean1 = 0;

  /* test ivariance */

  var1 = ivariance(groupone, 20);
  
  if (var1 > 14.42 && var1 < 14.425){
    printf("success: ivariance\n");
  }
  else{
    printf("ERROR: ivariance (calculated %f)\n", var1);
  }

  var1 = 0;

  /* test isd */

  sd1 = isd(groupone, 20);
  if (sd1 > 3.79 && sd1 < 3.78){
    printf("success: isd\n");
  }
  else{
    printf("ERROR: isd (calculated %f)\n", sd1);
  }

  sd1 = 0;

  /* test iipvariance */

  pv = iipvariance(groupone, grouptwo, 20, 20);
  if (pv > 18.84 && pv < 18.85 ){
    printf("success: iipvariance\n");
  }
  else{
    printf("ERROR: iipvariance (calculated %f)\n", pv);
  }

  pv = 0;

  /* test iittest */

  return 0;
}





