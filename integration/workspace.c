#include <config.h>
#include <stdlib.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

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

  w->order = (size_t *) malloc (n * sizeof (size_t));

  if (w->order == 0)
    {
      free (w->elist);
      free (w->rlist);
      free (w->blist);
      free (w->alist);
      free (w);		/* exception in constructor, avoid memory leak */

      GSL_ERROR_RETURN ("failed to allocate space for order ranges",
			GSL_ENOMEM, 0);
    }

  w->level = (size_t *) malloc (n * sizeof (size_t));

  if (w->level == 0)
    {
      free (w->order);
      free (w->elist);
      free (w->rlist);
      free (w->blist);
      free (w->alist);
      free (w);		/* exception in constructor, avoid memory leak */

      GSL_ERROR_RETURN ("failed to allocate space for order ranges",
			GSL_ENOMEM, 0);
    }

  w->size = 0 ;
  w->limit = n ;
  w->maximum_level = 0 ;
  
  return w ;
}

void
gsl_integration_workspace_free (gsl_integration_workspace * w)
{
  free (w->level) ;
  free (w->order) ;
  free (w->elist) ;
  free (w->rlist) ;
  free (w->blist) ;
  free (w->alist) ;
  free (w) ;
}

/*
size_t 
gsl_integration_workspace_limit (gsl_integration_workspace * w) 
{
  return w->limit ;
}


size_t 
gsl_integration_workspace_size (gsl_integration_workspace * w) 
{
  return w->size ;
}
*/
