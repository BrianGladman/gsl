#include <config.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_histogram.h>

int
gsl_histogram_increment (gsl_histogram * h, double x)
{
  int status = gsl_histogram_accumulate (h, x, 1.0);
  return status;
}

int
gsl_histogram_accumulate (gsl_histogram * h, double x, double weight)
{
  const size_t n = h->n;
  size_t index = 0;

  int status = gsl_histogram_find_impl (h->n, h->range, x, &index);

  if (status)
    {
      return GSL_EDOM;
    }

  if (index > n)
    {
      GSL_ERROR ("index lies outside valid range of 0 .. n - 1",
		 GSL_ESANITY);
    }

  h->bin[index] += weight;

  return GSL_SUCCESS;
}
