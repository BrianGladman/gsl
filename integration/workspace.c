#include <stdlib.h>
#include <gsl_integration.h>
#include <gsl_errno.h>

gsl_integration_workspace *
gsl_integration_workspace_alloc (const size_t n) 
{
  gsl_integration_workspace * w ;
  
  if (n == 0)
    {
      GSL_ERROR_RETURN ("workspace length n must be positive integer",
			GSL_EDOM, 0);
    }

  w = (gsl_integration_workspace *) 
    malloc (sizeof (gsl_integration_workspace));

  if (w == 0)
    {
      GSL_ERROR_RETURN ("failed to allocate space for workspace struct",
			GSL_ENOMEM, 0);
    }

  w->alist = (double *) malloc (n * sizeof (double));

  if (w->alist == 0)
    {
      free (w);		/* exception in constructor, avoid memory leak */

      GSL_ERROR_RETURN ("failed to allocate space for alist ranges",
			GSL_ENOMEM, 0);
    }

  w->blist = (double *) malloc (n * sizeof (double));

  if (w->blist == 0)
    {
      free (w->alist);
      free (w);		/* exception in constructor, avoid memory leak */

      GSL_ERROR_RETURN ("failed to allocate space for blist ranges",
			GSL_ENOMEM, 0);
    }

  w->rlist = (double *) malloc (n * sizeof (double));

  if (w->rlist == 0)
    {
      free (w->blist);
      free (w->alist);
      free (w);		/* exception in constructor, avoid memory leak */

      GSL_ERROR_RETURN ("failed to allocate space for rlist ranges",
			GSL_ENOMEM, 0);
    }


  w->elist = (double *) malloc (n * sizeof (double));

  if (w->elist == 0)
    {
      free (w->rlist);
      free (w->blist);
      free (w->alist);
      free (w);		/* exception in constructor, avoid memory leak */

      GSL_ERROR_RETURN ("failed to allocate space for elist ranges",
			GSL_ENOMEM, 0);
    }

  w->iord = (size_t *) malloc (n * sizeof (size_t));

  if (w->iord == 0)
    {
      free (w->elist);
      free (w->rlist);
      free (w->blist);
      free (w->alist);
      free (w);		/* exception in constructor, avoid memory leak */

      GSL_ERROR_RETURN ("failed to allocate space for iord ranges",
			GSL_ENOMEM, 0);
    }


  w->limit = n ;
  
  return w ;
}

void
gsl_integration_workspace_free (gsl_integration_workspace * w)
{
  free (w->iord) ;
  free (w->elist) ;
  free (w->rlist) ;
  free (w->blist) ;
  free (w->alist) ;
  free (w) ;
}
