#include "factorize.h"

int
FUNCTION(gsl_fft_halfcomplex,init) (size_t n,
				    TYPE(gsl_fft_wavetable_halfcomplex) * wavetable)
{
  int status;
  size_t n_factors;

  if (n == 0)
    {
      GSL_ERROR ("length n must be positive integer", GSL_EDOM);
    }

  wavetable->n = n;

  status = fft_halfcomplex_factorize (n, &n_factors, wavetable->factor);

  if (status)
    {
      GSL_ERROR ("factorization failed", GSL_EFACTOR);
    }

  wavetable->nf = n_factors;

  status = FUNCTION(gsl_fft_halfcomplex,generate) (n, wavetable);

  if (status)
    {
      GSL_ERROR ("could not generate wavetable", GSL_EFAILED);
    }

  return 0;
}

int
FUNCTION(gsl_fft_halfcomplex,generate) (size_t n,
						  TYPE(gsl_fft_wavetable_halfcomplex) * wavetable)
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
	      GSL_REAL(wavetable->trig[t]) = cos (theta);
	      GSL_IMAG(wavetable->trig[t]) = sin (theta);

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

TYPE(gsl_fft_wavetable_halfcomplex) *
FUNCTION(gsl_fft_halfcomplex,alloc) (size_t n)
{
  TYPE(gsl_fft_wavetable_halfcomplex) * w ;

  if (n == 0)
    {
      GSL_ERROR_RETURN ("length n must be positive integer", GSL_EDOM, 0);
    }

  w = (TYPE(gsl_fft_wavetable_halfcomplex) *) 
    malloc(sizeof(TYPE(gsl_fft_wavetable_halfcomplex)));

  if (w == NULL)
    {
      GSL_ERROR_RETURN ("failed to allocate struct", GSL_ENOMEM, 0);
    }

  w->scratch = (BASE *) malloc (2 * n * sizeof (BASE));

  if (w->scratch == NULL)
    {
      free(w) ; /* error in constructor, prevent memory leak */

      GSL_ERROR_RETURN ("failed to allocate scratch space", GSL_ENOMEM, 0);
    }

  w->trig = (gsl_complex *) malloc (n * sizeof (gsl_complex));

  if (w->trig == NULL)
    {
      free(w->scratch) ; /* error in constructor, prevent memory leak */
      free(w) ; 

      GSL_ERROR_RETURN ("failed to allocate trigonometric lookup table", 
			GSL_ENOMEM, 0);
    }

  return w;
}

void
FUNCTION(gsl_fft_halfcomplex,free) (TYPE(gsl_fft_wavetable_halfcomplex) * wavetable)
{

  /* release scratch space and trigonometric lookup tables */

  free (wavetable->scratch);
  wavetable->scratch = NULL;

  free (wavetable->trig);
  wavetable->trig = NULL;

  free (wavetable);
}
