#include <config.h>
#include <math.h>
#include <gsl_math.h>
#include <gsl_complex.h>
#include <gsl_poly.h>
#include <gsl_errno.h>

gsl_poly_complex_workspace * 
gsl_poly_complex_workspace_alloc (size_t n)
{
  size_t nc ;

  gsl_poly_complex_workspace * w ;
  
  if (n == 0)
    {
      GSL_ERROR_RETURN ("matrix size n must be positive integer", GSL_EDOM, 0);
    }

  w = (gsl_poly_complex_workspace *) 
    malloc (sizeof(gsl_poly_complex_workspace));

  if (w == 0)
    {
      GSL_ERROR_RETURN ("failed to allocate space for struct", GSL_ENOMEM, 0);
    }

  nc = n - 1;

  w->nc = nc;

  w->matrix = (double *) malloc (nc * nc * sizeof(double));

  if (w->matrix == 0)
    {
      free (w) ;       /* error in constructor, avoid memory leak */
      
      GSL_ERROR_RETURN ("failed to allocate space for workspace matrix", 
                        GSL_ENOMEM, 0);
    }

  return w ;
}

void 
gsl_poly_complex_workspace_free (gsl_poly_complex_workspace * w)
{
  free(w->matrix) ;
  free(w);
}

#include <stdio.h>

#define INDEX(m,i,j,n) ((m)[(i)*(n) + (j)])

void
dump (gsl_poly_complex_workspace * w)
{
  int i, j;

  for (i = 0; i < w->nc ; i++)
    {
      for (j = 0; j < w->nc ; j++)
        printf(" %10g", INDEX(w->matrix,i,j,w->nc)) ;
      
      printf("\n");
    }
  printf("\n");
}
