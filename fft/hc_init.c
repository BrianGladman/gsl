#include <config.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>

#include <gsl_errno.h>
#include <gsl_complex.h>
#include <gsl_fft_halfcomplex.h>

#include "factorize.h"

int
gsl_fft_halfcomplex_init (size_t n,
			  gsl_fft_halfcomplex_wavetable * wavetable)
{
  int status;
  size_t n_factors;

  if (n == 0)
    {
      GSL_ERROR ("length n must be positive integer", GSL_EDOM);
    }

  wavetable->n = n;

  status = gsl_fft_halfcomplex_factorize (n, &n_factors, wavetable->factor);

  if (status)
    {
      GSL_ERROR ("factorization failed", GSL_EFACTOR);
    }

  wavetable->nf = n_factors;

  status = gsl_fft_halfcomplex_generate_wavetable (n, wavetable);

  if (status)
    {
      GSL_ERROR ("could not generate wavetable", GSL_EFAILED);
    }

  return 0;
}

int
gsl_fft_halfcomplex_generate_wavetable (size_t n,
				  gsl_fft_halfcomplex_wavetable * wavetable)
{
  size_t i;
  double d_theta;
  size_t t, product, product_1, q;

  if (n == 0)
    {
      GSL_ERROR ("length n must be positive integer", GSL_EDOM);
    }

  if (n != wavetable->n)
    {
      GSL_ERROR ("wavetable does not match length of data", GSL_EINVAL);
    }

  d_theta = 2.0 * M_PI / ((double) n);

  t = 0;
  product = 1;
  for (i = 0; i < wavetable->nf; i++)
    {
      size_t j;
      const size_t factor = wavetable->factor[i];
      wavetable->twiddle[i] = wavetable->trig + t;
      product_1 = product;	/* product_1 = p_(i-1) */
      product *= factor;
      q = n / product;

      for (j = 1; j < factor; j++)
	{
	  size_t k;
	  size_t m = 0;
	  for (k = 1; k < (q + 1) / 2; k++)
	    {
	      double theta;
	      m = m + j * product_1;
	      m = m % n;
	      theta = d_theta * m;	/*  d_theta*j*k*product_1 */
	      wavetable->trig[t].real = cos (theta);
	      wavetable->trig[t].imag = sin (theta);

	      t++;
	    }
	}
    }

  if (t > (n / 2))
    {
      GSL_ERROR ("overflowed trigonometric lookup table", GSL_ESANITY);
    }

  return 0;
}

int
gsl_fft_halfcomplex_wavetable_alloc (size_t n,
				  gsl_fft_halfcomplex_wavetable * wavetable)
{
  if (n == 0)
    {
      GSL_ERROR ("length n must be positive integer", GSL_EDOM);
    }

  wavetable->scratch = malloc (n * sizeof (gsl_complex));

  if (wavetable->scratch == NULL)
    {
      GSL_ERROR ("failed to allocate scratch space", GSL_ENOMEM);
    }

  wavetable->trig = malloc (n * sizeof (gsl_complex));

  if (wavetable->trig == NULL)
    {
      GSL_ERROR ("failed to allocate trigonometric lookup table", GSL_ENOMEM);
    }

  return 0;
}

int
gsl_fft_halfcomplex_wavetable_free (gsl_fft_halfcomplex_wavetable * wavetable)
{

  /* release scratch space and trigonometric lookup tables */

  free (wavetable->scratch);
  wavetable->scratch = NULL;

  free (wavetable->trig);
  wavetable->trig = NULL;

  return 0;
}
