#include <config.h>
#include <stdlib.h>
#include <gsl_integration.h>
#include <gsl_errno.h>

gsl_integration_qawo_workspace *
gsl_integration_qawo_workspace_alloc (double omega, double L, 
				      enum gsl_integration_qawo_enum sine,
				      size_t n)
{
  gsl_integration_qawo_workspace *t;
  double * chebmo;

  t = (gsl_integration_qawo_workspace *)
    malloc (sizeof (gsl_integration_qawo_workspace));

  if (t == 0)
    {
      GSL_ERROR_RETURN ("failed to allocate space for qawo_workspace struct",
			GSL_ENOMEM, 0);
    }

  chebmo = (double *)  malloc (25 * n * sizeof (double));

  if (chebmo == 0)
    {
      free (t);
      GSL_ERROR_RETURN ("failed to allocate space for chebmo block",
			GSL_ENOMEM, 0);
    }

  t->i = 0;
  t->n = n;
  t->sine = sine;
  t->omega = omega;
  t->L = L;
  t->par = 0.5 * omega * L;
  t->chebmo = chebmo;

  return t;
}

int
gsl_integration_qawo_workspace_set (gsl_integration_qawo_workspace * t,
				    double omega, double L,
				    enum gsl_integration_qawo_enum sine)
{
  t->i = 0;
  t->omega = omega;
  t->sine = sine;
  t->L = L;
  t->par = 0.5 * omega * L;

  return GSL_SUCCESS;
}


int
gsl_integration_qawo_workspace_set_length (gsl_integration_qawo_workspace * t,
					   double L)
{
  /* return immediately if the length is the same as the old length */

  if (L == t->L)
    return GSL_SUCCESS;

  /* otherwise reset the table and compute the new parameters */

  t->i = 0;
  t->L = L;
  t->par = 0.5 * t->omega * L;

  return GSL_SUCCESS;
}


void
gsl_integration_qawo_workspace_free (gsl_integration_qawo_workspace * t)
{
  free (t->chebmo);
  free (t);
}

