#include <config.h>
#include <string.h>
#include <gsl_ieee_utils.h>
#include <gsl_test.h>

int 
main (void)
{
  {
    double d = 2.1 ; 
    const char mantissa[] 
      = "0000110011001100110011001100110011001100110011001101" ;
    gsl_ieee_double_rep r ;
    gsl_ieee_double_to_rep(&d, &r) ;

    gsl_test_int (r.sign,0,"x = 2.1, sign is +") ;
    gsl_test_int (r.exponent,1, "x = 2.1, exponent is 1") ;
    gsl_test_str (r.mantissa, mantissa,"x = 2.1, mantissa") ;
    gsl_test_int (r.type,GSL_IEEE_TYPE_NORMAL, "x = 2.1, type is NORMAL") ;
  }

  {
    double d = -1.3303577090924210146738460025517269968986511230468750 ; 
    const char mantissa[] 
      = "0101010010010010010100101010010010001000100011101110" ;
    gsl_ieee_double_rep r ;
    gsl_ieee_double_to_rep(&d, &r) ;

    gsl_test_int (r.sign, 1,"x = -1.3304..., sign is -") ;
    gsl_test_int (r.exponent,0, "x = -1.3304..., exponent is 0") ;
    gsl_test_str (r.mantissa, mantissa,"x = -1.3304..., mantissa") ;
    gsl_test_int (r.type,GSL_IEEE_TYPE_NORMAL, 
		  "x = -1.3304..., type is NORMAL") ;
  }

  {
    double d = 3.37e297 ;
    const char mantissa[] 
      = "0100100111001001100101111001100000100110011101000100" ;
    gsl_ieee_double_rep r ;
    gsl_ieee_double_to_rep(&d, &r) ;

    gsl_test_int (r.sign, 0,"x = 3.37e297..., sign is +") ;
    gsl_test_int (r.exponent, 988, "x = 3.37e297..., exponent is 998") ;
    gsl_test_str (r.mantissa, mantissa,"x = 3.37e297..., mantissa") ;
    gsl_test_int (r.type,GSL_IEEE_TYPE_NORMAL, 
		  "x = 3.37e297..., type is NORMAL") ;
  }

  return gsl_test_summary() ;
}

