#include <config.h>
#include <gsl_errno.h>
#include <gsl_histogram.h>

void
gsl_histogram_reset (gsl_histogram * h)
{
  size_t i;
  const size_t n = h->n;

  for (i = 0; i < n; i++)
    {
      h->bin[i] = 0;
    }
}
