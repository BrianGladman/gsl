double
FUNCTION(gsl_stats,wabsdev) (const BASE w[], const size_t wstride, const BASE data[], const size_t stride, const size_t n)
{
  const double wmean = FUNCTION(gsl_stats,wmean)(w, wstride, data, stride, n);
  return FUNCTION(gsl_stats,wabsdev_m)(w, wstride, data, stride, n, wmean);
}
    
double 
FUNCTION(gsl_stats,wabsdev_m) (const BASE w[], const size_t wstride, const BASE data[], const size_t stride, const size_t n, const double wmean)
{
  /* Compute the weighted absolute deviation of a dataset */

  long double wabsdev = 0;
  long double W = 0;

  size_t i;

  /* find the sum of the absolute deviations */
  for (i = 0; i < n; i++)
    {
      BASE wi = w[i * wstride];
      
      if (wi > 0) {
        const long double delta = fabs(data[i * stride] - wmean);
        W += wi ;
        wabsdev += (delta - wabsdev) * (wi / W);
      }
    }

  return wabsdev;
}

