#include <config.h>
#include <gsl_errno.h>
#include <gsl_histogram.h>

double
gsl_histogram_max (const gsl_histogram * h)
{
  const int n = h->n;

  return h->range[n];
}

double
gsl_histogram_min (const gsl_histogram * h)
{
  return h->range[0];
}

size_t
gsl_histogram_bins (const gsl_histogram * h)
{
  return h->n;
}
