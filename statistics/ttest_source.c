double 
FUNCTION(gsl_stats,ttest) (const BASE data1[], 
                           const size_t stride1, const size_t n1, 
                           const BASE data2[],
                           const size_t stride2, const size_t n2)
{
  /* runs a t-test between two datasets representing independent
     samples. Tests to see if the difference between means of the
     samples is different from zero */

  /* find means for the two samples */
  const double mean1 = FUNCTION(gsl_stats,mean) (data1, stride1, n1);
  const double mean2 = FUNCTION(gsl_stats,mean) (data2, stride2, n2);

  /* find pooled variance for the two samples */
  const double pv = FUNCTION(gsl_stats,pvariance) (data1, stride1, n1, data2, stride2, n2);

  /* calculate the t statistic */
  const double t = (mean1 - mean2) / (sqrt (pv * ((1.0 / n1) + (1.0 / n2))));

  return t;
}

