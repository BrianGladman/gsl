#include <config.h>
#include <stdlib.h>
#include <string.h>
#include <gsl_ieee_utils.h>
#include <gsl_errno.h>

static int 
lookup_string (const char * p, int * precision, int * rounding, 
	       int * exception_mask) ;


void
gsl_ieee_env_setup (void)
{
  const char * p = getenv("GSL_IEEE_MODE") ;
  
  if (p)
    gsl_ieee_read_mode_string (p) ;

}

int
gsl_ieee_read_mode_string (const char * description)
{
  char * start = strdup (description) ;
  char * end;

  char * p = start;

  int precision = 0 ;
  int precision_count = 0 ;

  int rounding = 0 ;
  int rounding_count = 0 ;

  int exception_mask = 0 ;
  int exception_count = 0 ;

  do {
    int status ;
    int new_precision, new_rounding, new_exception ;

    end = strchr (p,';') ;

    if (end) 
      {
	*end = '\0' ;
	do 
	  {
	    end++ ;  /* skip over trailing whitespace */
	  } 
	while (*end == ' ' || *end == ';') ;
      }
	
    new_precision = 0 ; 
    new_rounding = 0 ; 
    new_exception = 0 ;

    status = lookup_string (p, &new_precision, &new_rounding, &new_exception) ;

    printf("%s: %d, %d, %d\n", p, new_precision,new_rounding,new_exception) ;

    if (status)
      GSL_ERROR ("unrecognized IEEE mode. Valid settings are:\n\n" 
		 "  single double extended\n"
		 "  round-to-nearest round-down round-up round-to-zero\n"
		 "  mask-invalid mask-denormalized mask-overflow "
		 "mask-underflow catch-inexact\n"
		 "\nseparated by semicolons. "
		 "(default is GSL_IEEE_MODE=\"extended;round-to-nearest\")",
		 GSL_EINVAL) ;

    if (new_precision) 
      {
	precision = new_precision ;
	precision_count ++ ;
	if (precision_count > 1)
	  GSL_ERROR ("attempted to set IEEE precision twice", GSL_EINVAL) ;
      }

    if (new_rounding) 
      {
	rounding = new_rounding ;
	rounding_count ++ ;
	if (rounding_count > 1)
	  GSL_ERROR ("attempted to set IEEE rounding mode twice", GSL_EINVAL) ;
      }

    if (new_exception) 
      {
	exception_mask |= new_exception ;
	exception_count ++ ;
      }

    p = end ; 

  } while (end && *p != '\0') ;

  gsl_ieee_set_mode (precision, rounding, exception_mask) ;

  free(start) ;

  return GSL_SUCCESS ;
}

static int 
lookup_string (const char * p, int * precision, int * rounding, 
	       int * exception_mask)
{
  if (strcmp(p,"single") == 0) 
    {
      *precision = GSL_IEEE_SINGLE_PRECISION ;
    }
  else if (strcmp(p,"double") == 0) 
    {
      *precision = GSL_IEEE_DOUBLE_PRECISION ;
    }
  else if (strcmp(p,"extended") == 0) 
    {
      *precision = GSL_IEEE_EXTENDED_PRECISION ;
    }
  else if (strcmp(p,"round-to-nearest") == 0) 
    {
      *rounding = GSL_IEEE_ROUND_TO_NEAREST ;
    }
  else if (strcmp(p,"round-down") == 0) 
    {
      *rounding = GSL_IEEE_ROUND_DOWN ;
    }
  else if (strcmp(p,"round-up") == 0) 
    {
      *rounding = GSL_IEEE_ROUND_UP ;
    }
  else if (strcmp(p,"round-to-zero") == 0) 
    {
      *rounding = GSL_IEEE_ROUND_TO_ZERO ;
    }
  else if (strcmp(p,"mask-invalid") == 0) 
    {
      *exception_mask = GSL_IEEE_MASK_INVALID ;
    }
  else if (strcmp(p,"mask-denormalized") == 0) 
    {
      *exception_mask = GSL_IEEE_MASK_DENORMALIZED ;
    }
  else if (strcmp(p,"mask-division-by-zero") == 0) 
    {
      *exception_mask = GSL_IEEE_MASK_DIVISION_BY_ZERO ;
    }
  else if (strcmp(p,"mask-overflow") == 0) 
    {
      *exception_mask = GSL_IEEE_MASK_OVERFLOW ;
    }
  else if (strcmp(p,"mask-underflow") == 0) 
    {
      *exception_mask = GSL_IEEE_MASK_UNDERFLOW ;
    }
  else if (strcmp(p,"catch-inexact") == 0) 
    {
      *exception_mask = GSL_IEEE_CATCH_INEXACT ;
    }
  else 
    {
      return 1 ;
    }

  return 0 ;
}




