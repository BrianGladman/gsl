#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_legendre.h>

int
main(void)
{
  const size_t lmax = 5;
  const size_t plm_size = gsl_sf_legendre_array_n(lmax);
  const size_t nlm = gsl_sf_legendre_nlm(lmax);
  const size_t idx21 = gsl_sf_legendre_array_index(2, 1);
  double * Plm = malloc(plm_size * sizeof(double));
  double * dPlm = malloc(nlm * sizeof(double));
  double * d2Plm = malloc(nlm * sizeof(double));
  double x;

  gsl_sf_legendre_precompute(GSL_SF_LEGENDRE_SPHARM, lmax, 1, Plm);

  for (x = -1.0; x <= 1.0; x += 0.01)
    {
      gsl_sf_legendre_deriv2_alt_arrayx(GSL_SF_LEGENDRE_SPHARM, lmax, x, Plm, dPlm, d2Plm);
      printf("%f %e %e %e\n", x, Plm[idx21], dPlm[idx21], d2Plm[idx21]);
    }

  free(Plm);
  free(dPlm);
  free(d2Plm);

  return 0;
}
