double
FUNCTION (gsl_stats, wmean) (const BASE w[], size_t wstride, const BASE data[], const size_t stride, const size_t size)
{
  /* Compute the weighted arithmetic mean M of a dataset using the
     recurrence relation

     M(n) = M(n-1) + (data[n] - M(n-1)) (w(n)/(W(n-1) + w(n))) 
     W(n) = W(n-1) + w(n)

   */

  long double wmean = 0;
  long double W = 0;

  size_t i;

  for (i = 0; i < size; i++)
    {
      BASE wi = w[i * wstride];

      if (wi > 0)
        {
          W += wi;
          wmean += (data[i * stride] - wmean) * (wi / W);
        }
    }

  return wmean;
}
