double 
FUNCTION(gsl_stats,wskew) (const BASE w[], const size_t wstride, const BASE data[], const size_t stride, const size_t n)
{
  const double wmean = FUNCTION(gsl_stats,wmean)(w, wstride, data, stride, n);
  const double wsd = FUNCTION(gsl_stats,wsd_m)(w, wstride, data, stride, n, wmean);
  return FUNCTION(gsl_stats,wskew_m_sd)(w, wstride, data, stride, n, wmean, wsd);
}
    
double 
FUNCTION(gsl_stats,wskew_m_sd) (const BASE w[], const size_t wstride, 
                                const BASE data[], 
                                const size_t stride, const size_t n,
                                const double wmean, const double wsd)
{
  /* Compute the weighted skewness of a dataset */

  long double wskew = 0;
  long double W = 0;

  size_t i;

  /* find the sum of the cubed deviations, normalized by the sd. */

  /* we use a recurrence relation to stably update a running value so
     there aren't any large sums that can overflow */

  for (i = 0; i < n; i++)
    {
      BASE wi = w[i * wstride];
      
      if (wi > 0) {
        const long double x = (data[i * stride] - wmean) / wsd;
        W += wi ;
        wskew += (x * x * x - wskew) * (wi / W);
      }
    }

  return wskew;
}

