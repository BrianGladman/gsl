#include <config.h>
#include <stdlib.h>
#include <gsl_integration.h>
#include <gsl_errno.h>

gsl_integration_qawf_workspace *
gsl_integration_qawf_workspace_alloc (double par, size_t n)
{
  gsl_integration_qawf_workspace *t;
  double * chebmo;

  t = (gsl_integration_qawf_workspace *)
    malloc (sizeof (gsl_integration_qawf_workspace));

  if (t == 0)
    {
      GSL_ERROR_RETURN ("failed to allocate space for qawf_workspace struct",
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
  t->par = par;
  t->chebmo = chebmo;

  return t;
}

int
gsl_integration_qawf_workspace_set (gsl_integration_qawf_workspace * t,
				double par)
{
  t->i = 0;
  t->par = par;

  return GSL_SUCCESS;
}


void
gsl_integration_qawf_workspace_free (gsl_integration_qawf_workspace * t)
{
  free (t->chebmo);
  free (t);
}

