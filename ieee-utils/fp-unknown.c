#include <config.h>
#include <gsl_ieee_utils.h>

void
gsl_ieee_set_mode (int precision, int rounding, int exception_mask)
{
  
  GSL_ERROR ("the IEEE interface for this platform could not be determined "
	     "at configure time\n",
	     GSL_EUNSUP)
  
}
