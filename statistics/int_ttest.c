/* This is an automatically generated file created by
   convert_double_to_int.pl -- do not edit                */

#include <math.h>
#include <gsl_statistics_int.h>

double 
gsl_stats_int_ttest (const int data1[], const int data2[],
		 const size_t n1, const size_t n2)
{
  /* runs a t-test between two integer datasets representing independent
     samples. Tests to see if the difference between means of the
     samples is different from zero */

  /* find means for the two samples */
  const double mean1 = gsl_stats_int_mean (data1, n1);
  const double mean2 = gsl_stats_int_mean (data2, n2);

  /* find pooled variance for the two samples */
  const double pv = gsl_stats_int_pvariance (data1, data2, n1, n2);

  /* calculate the t statistic */
  const double t = (mean1 - mean2) / (sqrt (pv * ((1.0 / n1) + (1.0 / n2))));

  return t;
}

