static double 
FUNCTION(compute,wvariance) (const BASE w[], const size_t wstride, const BASE data[], const size_t stride, const size_t n, const double wmean);

static double
FUNCTION(compute,factor) (const BASE w[], const size_t wstride, const size_t n);

static double
FUNCTION(compute,wvariance) (const BASE w[], const size_t wstride, const BASE data[], const size_t stride, const size_t n, const double wmean)
{
  /* takes a dataset and finds the weighted variance */

  long double wvariance = 0 ;
  long double W = 0;

  size_t i;

  /* find the sum of the squares */
  for (i = 0; i < n; i++)
    {
      BASE wi = w[i * wstride];

      if (wi > 0) {
        const long double delta = (data[i * stride] - wmean);
        W += wi ;
        wvariance += (delta * delta - wvariance) * (wi / W);
      }
    }

  return wvariance ;
}

static double
FUNCTION(compute,factor) (const BASE w[], const size_t wstride, const size_t n)
{
  /* Find the factor ``N/(N-1)'' which multiplies the raw std dev */

  long double a = 0 ;
  long double b = 0;
  long double factor;

  size_t i;

  /* find the sum of the squares */
  for (i = 0; i < n; i++)
    {
      BASE wi = w[i * wstride];

      if (wi > 0)
        {
          a += wi ;
          b += wi * wi ;
        }
    }

  factor = (a*a) / ((a*a) - b);

  return factor ;
}

double 
FUNCTION(gsl_stats,wvariance_with_mean) (const BASE w[], const size_t wstride, const BASE data[], const size_t stride, const size_t n, const double wmean)
{
  const double wvariance = FUNCTION(compute,wvariance) (w, wstride, data, stride, n, wmean);
  return wvariance;
}

double 
FUNCTION(gsl_stats,est_wvariance_with_mean) (const BASE w[], const size_t wstride, const BASE data[], const size_t stride, const size_t n, const double wmean)
{
  const double variance = FUNCTION(compute,wvariance) (w, wstride, data, stride, n, wmean);
  const double scale = FUNCTION(compute,factor)(w, wstride, n);
  const double est_variance = scale * variance ;
  
  return est_variance;
}

double 
FUNCTION(gsl_stats,wsd_with_mean) (const BASE w[], const size_t wstride, const BASE data[], const size_t stride, const size_t n, const double wmean)
{
  const double wvariance = FUNCTION(compute,wvariance) (w, wstride, data, stride, n, wmean);
  const double wsd = sqrt (wvariance);

  return wsd;
}

double 
FUNCTION(gsl_stats,est_wsd_with_mean) (const BASE w[], const size_t wstride, const BASE data[], const size_t stride, const size_t n, const double wmean)
{
  const double variance = FUNCTION(compute,wvariance) (w, wstride, data, stride, n, wmean);
  const double scale = FUNCTION(compute,factor)(w, wstride, n);
  const double est_wsd = sqrt(scale * variance) ;
  
  return est_wsd;
}

double 
FUNCTION(gsl_stats,wvariance) (const BASE w[], const size_t wstride, const BASE data[], const size_t stride, const size_t n)
{
  const double wmean = FUNCTION(gsl_stats,wmean) (w, wstride, data, stride, n);
  return FUNCTION(gsl_stats,wvariance_with_mean)(w, wstride, data, stride, n, wmean);
}

double 
FUNCTION(gsl_stats,wsd) (const BASE w[], const size_t wstride, const BASE data[], const size_t stride, const size_t n)
{
  const double wmean = FUNCTION(gsl_stats,wmean) (w, wstride, data, stride, n);
  return FUNCTION(gsl_stats,wsd_with_mean) (w, wstride, data, stride, n, wmean) ;
}

double 
FUNCTION(gsl_stats,est_wvariance) (const BASE w[], const size_t wstride, const BASE data[], const size_t stride, const size_t n)
{
  const double wmean = FUNCTION(gsl_stats,wmean) (w, wstride, data, stride, n);
  return FUNCTION(gsl_stats,est_wvariance_with_mean)(w, wstride, data, stride, n, wmean);
}

double 
FUNCTION(gsl_stats,est_wsd) (const BASE w[], const size_t wstride, const BASE data[], const size_t stride, const size_t n)
{
  const double wmean = FUNCTION(gsl_stats,wmean) (w, wstride, data, stride, n);
  return FUNCTION(gsl_stats,est_wsd_with_mean) (w, wstride, data, stride, n, wmean) ;
}

