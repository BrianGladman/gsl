#include <math.h>
#include <gsl_statistics.h>

double 
gsl_stats_ittest (const int data1[], const int data2[],
		 unsigned int n1, unsigned int n2)
{
  /* runs a t-test between two integer datasets representing independent
     samples. Tests to see if the difference between means of the
     samples is different from zero */

  double mean1, mean2;		/* means of the two samples */
  double sd1, sd2;		/* standard deviations */
  double pv;			/* pooled variance */
  double t;			/* the t statistic */

  /* find means and standard deviations for the two samples */
  mean1 = gsl_stats_imean (data1, n1);
  mean2 = gsl_stats_imean (data2, n2);
  sd1 = gsl_stats_iest_stddev (data1, n1);
  sd2 = gsl_stats_iest_stddev (data2, n2);
  pv = gsl_stats_ipvariance (data1, data2, n1, n2);

  /* calculate the t statistic */
  t = (mean1 - mean2) / (sqrt (pv * ((1.0 / n1) + (1.0 / n2))));

  return t;
}

