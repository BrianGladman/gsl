#include <math.h>
#include <float.h>
#include <gsl_integration.h>

#include "err.h"

double rescale_error (double err, 
		      const double result_abs, const double result_asc)
{
  err = fabs(err) ;

  if (result_asc != 0 && err != 0)
      {
	double scale = pow((200 * err / result_asc), 1.5) ;
	
	if (scale < 1)
	  {
	    err = result_asc * scale ;
	  }
	else 
	  {
	    err = result_asc ;
	  }
      }
  if (result_abs > DBL_MIN / (50 * DBL_EPSILON))
    {
      double min_err = (DBL_EPSILON * 50) * result_abs ;

      if (min_err > err) 
	{
	  err = min_err ;
	}
    }
  
  return err ;
}

