#include <config.h>
#include <stdlib.h>
#include <gsl_integration.h>
#include <gsl_errno.h>

gsl_integration_workspace_pts *
gsl_integration_workspace_pts_alloc (const size_t npts) 
{
  gsl_integration_workspace_pts * w ;
  
  if (npts == 0)
    {
      GSL_ERROR_RETURN ("workspace_pts length npts must be positive integer",
			GSL_EDOM, 0);
    }

  w = (gsl_integration_workspace_pts *) 
    malloc (sizeof (gsl_integration_workspace_pts));

  if (w == 0)
    {
      GSL_ERROR_RETURN ("failed to allocate space for workspace_pts struct",
			GSL_ENOMEM, 0);
    }

  w->level = (unsigned int *) malloc (npts * sizeof (unsigned int));

  if (w->level == 0)
    {
      free (w);		/* exception in constructor, avoid memory leak */

      GSL_ERROR_RETURN ("failed to allocate space for levels",
			GSL_ENOMEM, 0);
    }

  w->ndin = (unsigned int *) malloc (npts * sizeof (unsigned int));

  if (w->ndin == 0)
    {
      free (w->level);
      free (w);		/* exception in constructor, avoid memory leak */

      GSL_ERROR_RETURN ("failed to allocate space for ndin flags",
			GSL_ENOMEM, 0);
    }

  w->npts = npts ;

  return w ;
}

void
gsl_integration_workspace_pts_free (gsl_integration_workspace_pts * w)
{
  free (w->ndin) ;
  free (w->level) ;
  free (w) ;
}
