double 
FUNCTION(gsl_stats,wkurtosis) (const BASE w[], const size_t wstride, const BASE data[], const size_t stride, const size_t n)
{
  const double wmean = FUNCTION(gsl_stats,wmean)(w, wstride, data, stride, n);
  const double est_wsd = FUNCTION(gsl_stats,est_wsd_with_mean)(w, wstride, data, stride, n, wmean);
  return FUNCTION(gsl_stats,wkurtosis_with_mean_and_sd)(w, wstride, data, stride, n, wmean, est_wsd);
}
    
double 
FUNCTION(gsl_stats,wkurtosis_with_mean_and_sd) (const BASE w[], const size_t wstride,
                                                const BASE data[], 
                                                const size_t stride,
                                                const size_t n,
                                                const double wmean, 
                                                const double wsd)
{
  /* takes a dataset and finds the kurtosis */

  long double wavg = 0, kurtosis;
  long double W = 0;
  size_t i;

  /* find the fourth moment the deviations, normalized by the sd */

  /* we use a recurrence relation to stably update a running value so
     there aren't any large sums that can overflow */

  for (i = 0; i < n; i++)
    {
      BASE wi = w[i * wstride];
      
      if (wi > 0) {
        const long double x = (data[i * stride] - wmean) / wsd;
        W += wi ;
        wavg += (x * x * x * x - wavg) * (wi / W);
      }
    }

  kurtosis = wavg - 3.0;  /* makes kurtosis zero for a Gaussian */

  return kurtosis;
}

