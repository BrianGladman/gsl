#include <config.h>
#include <math.h>
#include <string.h>
#include <gsl_ieee_utils.h>
#include <gsl_test.h>

int
main (void)
{

  /* Check for +ZERO (float) */

  {
    float f = 0.0;
    const char mantissa[] = "00000000000000000000000";
    gsl_ieee_float_rep r;
    gsl_ieee_float_to_rep (&f, &r);

    gsl_test_int (r.sign, 0, "x = 0, sign is +");
    gsl_test_int (r.exponent, -127, "x = 0, exponent is -127");
    gsl_test_str (r.mantissa, mantissa, "x = 0, mantissa");
    gsl_test_int (r.type, GSL_IEEE_TYPE_ZERO, "x = 0, type is ZERO");
  }

  /* Check for -ZERO (float) */

  {
    float f = 0.0;
    const char mantissa[] = "00000000000000000000000";
    gsl_ieee_float_rep r;

    float x = f * -1.0;
    gsl_ieee_float_to_rep (&x, &r);

    gsl_test_int (r.sign, 1, "x = -1*0, sign is -");
    gsl_test_int (r.exponent, -127, "x = -1*0, exponent is -127");
    gsl_test_str (r.mantissa, mantissa, "x = -1*0, mantissa");
    gsl_test_int (r.type, GSL_IEEE_TYPE_ZERO, "x = -1*0, type is ZERO");
  }

  /* Check for a positive NORMAL number (e.g. 2.1) (float) */

  {
    float f = 2.1;
    const char mantissa[] = "00001100110011001100110";

    gsl_ieee_float_rep r;
    gsl_ieee_float_to_rep (&f, &r);

    gsl_test_int (r.sign, 0, "x = 2.1, sign is +");
    gsl_test_int (r.exponent, 1, "x = 2.1, exponent is 1");
    gsl_test_str (r.mantissa, mantissa, "x = 2.1, mantissa");
    gsl_test_int (r.type, GSL_IEEE_TYPE_NORMAL, "x = 2.1, type is NORMAL");
  }


  /* Check for a negative NORMAL number (e.g. -1.3304...) (float) */

  {
    float f = -1.3303577090924210 ;
    const char mantissa[] = "01010100100100100101001";

    gsl_ieee_float_rep r;
    gsl_ieee_float_to_rep (&f, &r);

    gsl_test_int (r.sign, 1, "x = -1.3304..., sign is -");
    gsl_test_int (r.exponent, 0, "x = -1.3304..., exponent is 0");
    gsl_test_str (r.mantissa, mantissa, "x = -1.3304..., mantissa");
    gsl_test_int (r.type, GSL_IEEE_TYPE_NORMAL,
		  "x = -1.3304..., type is NORMAL");
  }

  /* Check for a large positive NORMAL number (e.g. 3.37e31) (float) */

  {
    float f = 3.37e31;
    const char mantissa[] = "10101001010110101001001";
    gsl_ieee_float_rep r;
    gsl_ieee_float_to_rep (&f, &r);

    gsl_test_int (r.sign, 0, "x = 3.37e31, sign is +");
    gsl_test_int (r.exponent, 104, "x = 3.37e31, exponent is 104");
    gsl_test_str (r.mantissa, mantissa, "x = 3.37e31, mantissa");
    gsl_test_int (r.type, GSL_IEEE_TYPE_NORMAL, "x = 3.37e31, type is NORMAL");
  }

  /* Check for a small positive NORMAL number (e.g. 3.37e-31) (float) */

  {
    float f = 3.37e-31;
    const char mantissa[] = "10110101011100110111011";

    gsl_ieee_float_rep r;
    gsl_ieee_float_to_rep (&f, &r);

    gsl_test_int (r.sign, 0, "x = 3.37e-31, sign is +");
    gsl_test_int (r.exponent, -102, "x = 3.37e-31, exponent is -102");
    gsl_test_str (r.mantissa, mantissa, "x = 3.37e-31, mantissa");
    gsl_test_int (r.type, GSL_IEEE_TYPE_NORMAL,
		  "x = 3.37e-31, type is NORMAL");
  }

  /* Check for FLT_MIN (smallest possible number that is not denormal) */

  {
    float f = 1.17549435e-38;	/* FLT_MIN (float) */
    const char mantissa[] = "00000000000000000000000";
    gsl_ieee_float_rep r;
    gsl_ieee_float_to_rep (&f, &r);

    gsl_test_int (r.sign, 0, "x = FLT_MIN, sign is +");
    gsl_test_int (r.exponent, -126, "x = FLT_MIN, exponent is -126");
    gsl_test_str (r.mantissa, mantissa, "x = FLT_MIN, mantissa");
    gsl_test_int (r.type, GSL_IEEE_TYPE_NORMAL, "x = FLT_MIN, type is NORMAL");
  }

  /* Check for FLT_MAX (largest possible number that is not Inf) */

  {
    float f = 3.40282347e+38;	/* FLT_MAX */
    const char mantissa[] = "11111111111111111111111";

    gsl_ieee_float_rep r;
    gsl_ieee_float_to_rep (&f, &r);

    gsl_test_int (r.sign, 0, "x = FLT_MAX, sign is +");
    gsl_test_int (r.exponent, 127, "x = FLT_MAX, exponent is 127");
    gsl_test_str (r.mantissa, mantissa, "x = FLT_MAX, mantissa");
    gsl_test_int (r.type, GSL_IEEE_TYPE_NORMAL, "x = FLT_MAX, type is NORMAL");
  }

  /* Check for DENORMAL numbers (e.g. FLT_MIN/2^n) */

  {
    float f = 1.17549435e-38;	/* FLT_MIN */
    char mantissa[] = "10000000000000000000000";

    int i;
    gsl_ieee_float_rep r;

    for (i = 0; i < 23; i++)
      {
	float x = f / pow (2.0, 1 + (float) i);
	mantissa[i] = '1';
	gsl_ieee_float_to_rep (&x, &r);

	gsl_test_int (r.sign, 0, "x = FLT_MIN/2^%d, sign is +", i + 1);
	gsl_test_int (r.exponent, -127,
		      "x = FLT_MIN/2^%d, exponent is -127", i + 1);
	gsl_test_str (r.mantissa, mantissa,
		      "x = FLT_MIN/2^%d, mantissa", i + 1);
	gsl_test_int (r.type, GSL_IEEE_TYPE_DENORMAL,
		      "x = FLT_MIN/2^%d, type is DENORMAL", i + 1);
	mantissa[i] = '0';
      }
  }

  /* Check for positive INFINITY (e.g. 2*FLT_MAX) */

  {
    float f = 3.40282347e+38;	/* FLT_MAX */
    const char mantissa[] = "00000000000000000000000";

    gsl_ieee_float_rep r;

    float x = 2 * f;
    gsl_ieee_float_to_rep (&x, &r);

    gsl_test_int (r.sign, 0, "x = 2*FLT_MAX, sign is +");
    gsl_test_int (r.exponent, 128, "x = 2*FLT_MAX, exponent is 128");
    gsl_test_str (r.mantissa, mantissa, "x = 2*FLT_MAX, mantissa");
    gsl_test_int (r.type, GSL_IEEE_TYPE_INF, "x = -2*FLT_MAX, type is INF");
  }

  /* Check for negative INFINITY (e.g. -2*FLT_MAX) */

  {
    float f = 3.40282347e+38;	/* FLT_MAX */
    const char mantissa[] = "00000000000000000000000";

    gsl_ieee_float_rep r;

    float x = -2 * f;
    gsl_ieee_float_to_rep (&x, &r);

    gsl_test_int (r.sign, 1, "x = -2*FLT_MAX, sign is -");
    gsl_test_int (r.exponent, 128, "x = -2*FLT_MAX, exponent is 128");
    gsl_test_str (r.mantissa, mantissa, "x = -2*FLT_MAX, mantissa");
    gsl_test_int (r.type, GSL_IEEE_TYPE_INF, "x = -2*FLT_MAX, type is INF");
  }

  /* Check for NAN (e.g. Inf - Inf) (float) */

  {
    const char mantissa[] = "00000000000000000000000";
    gsl_ieee_float_rep r;
    float x, y, z;

    x = 1.0 / 0.0;
    y = 2.0 / 0.0;
    z = y - x;

    gsl_ieee_float_to_rep (&z, &r);

    /* We don't check the sign and we don't check the mantissa because
       they could be anything for a NaN */

    gsl_test_int (r.exponent, 128, "x = NaN, exponent is 128");
    gsl_test_int (r.type, GSL_IEEE_TYPE_NAN, "x = NaN, type is NAN");
  }


  /* Check for +ZERO */

  {
    double d = 0.0;
    const char mantissa[]
    = "0000000000000000000000000000000000000000000000000000";
    gsl_ieee_double_rep r;
    gsl_ieee_double_to_rep (&d, &r);

    gsl_test_int (r.sign, 0, "x = 0, sign is +");
    gsl_test_int (r.exponent, -1023, "x = 0, exponent is -1023");
    gsl_test_str (r.mantissa, mantissa, "x = 0, mantissa");
    gsl_test_int (r.type, GSL_IEEE_TYPE_ZERO, "x = 0, type is ZERO");
  }

  /* Check for -ZERO */

  {
    double d = 0.0;
    const char mantissa[]
    = "0000000000000000000000000000000000000000000000000000";
    gsl_ieee_double_rep r;

    double x = d * -1.0;
    gsl_ieee_double_to_rep (&x, &r);

    gsl_test_int (r.sign, 1, "x = -1*0, sign is -");
    gsl_test_int (r.exponent, -1023, "x = -1*0, exponent is -1023");
    gsl_test_str (r.mantissa, mantissa, "x = -1*0, mantissa");
    gsl_test_int (r.type, GSL_IEEE_TYPE_ZERO, "x = -1*0, type is ZERO");
  }

  /* Check for a positive NORMAL number (e.g. 2.1) */

  {
    double d = 2.1;
    const char mantissa[]
    = "0000110011001100110011001100110011001100110011001101";
    gsl_ieee_double_rep r;
    gsl_ieee_double_to_rep (&d, &r);

    gsl_test_int (r.sign, 0, "x = 2.1, sign is +");
    gsl_test_int (r.exponent, 1, "x = 2.1, exponent is 1");
    gsl_test_str (r.mantissa, mantissa, "x = 2.1, mantissa");
    gsl_test_int (r.type, GSL_IEEE_TYPE_NORMAL, "x = 2.1, type is NORMAL");
  }


  /* Check for a negative NORMAL number (e.g. -1.3304...) */

  {
    double d = -1.3303577090924210146738460025517269968986511230468750;
    const char mantissa[]
    = "0101010010010010010100101010010010001000100011101110";
    gsl_ieee_double_rep r;
    gsl_ieee_double_to_rep (&d, &r);

    gsl_test_int (r.sign, 1, "x = -1.3304..., sign is -");
    gsl_test_int (r.exponent, 0, "x = -1.3304..., exponent is 0");
    gsl_test_str (r.mantissa, mantissa, "x = -1.3304..., mantissa");
    gsl_test_int (r.type, GSL_IEEE_TYPE_NORMAL,
		  "x = -1.3304..., type is NORMAL");
  }

  /* Check for a large positive NORMAL number (e.g. 3.37e297) */

  {
    double d = 3.37e297;
    const char mantissa[]
    = "0100100111001001100101111001100000100110011101000100";
    gsl_ieee_double_rep r;
    gsl_ieee_double_to_rep (&d, &r);

    gsl_test_int (r.sign, 0, "x = 3.37e297, sign is +");
    gsl_test_int (r.exponent, 988, "x = 3.37e297, exponent is 998");
    gsl_test_str (r.mantissa, mantissa, "x = 3.37e297, mantissa");
    gsl_test_int (r.type, GSL_IEEE_TYPE_NORMAL,
		  "x = 3.37e297, type is NORMAL");
  }

  /* Check for a small positive NORMAL number (e.g. 3.37e-297) */

  {
    double d = 3.37e-297;
    const char mantissa[]
    = "0001101000011011101011100001110010100001001100110111";
    gsl_ieee_double_rep r;
    gsl_ieee_double_to_rep (&d, &r);

    gsl_test_int (r.sign, 0, "x = 3.37e-297, sign is +");
    gsl_test_int (r.exponent, -985, "x = 3.37e-297, exponent is -985");
    gsl_test_str (r.mantissa, mantissa, "x = 3.37e-297, mantissa");
    gsl_test_int (r.type, GSL_IEEE_TYPE_NORMAL,
		  "x = 3.37e-297, type is NORMAL");
  }

  /* Check for DBL_MIN (smallest possible number that is not denormal) */

  {
    double d = 2.2250738585072014e-308;		/* DBL_MIN */
    const char mantissa[]
    = "0000000000000000000000000000000000000000000000000000";
    gsl_ieee_double_rep r;
    gsl_ieee_double_to_rep (&d, &r);

    gsl_test_int (r.sign, 0, "x = DBL_MIN, sign is +");
    gsl_test_int (r.exponent, -1022, "x = DBL_MIN, exponent is -1022");
    gsl_test_str (r.mantissa, mantissa, "x = DBL_MIN, mantissa");
    gsl_test_int (r.type, GSL_IEEE_TYPE_NORMAL,
		  "x = DBL_MIN, type is NORMAL");
  }

  /* Check for DBL_MAX (largest possible number that is not Inf) */

  {
    double d = 1.7976931348623157e+308;		/* DBL_MAX */
    const char mantissa[]
    = "1111111111111111111111111111111111111111111111111111";
    gsl_ieee_double_rep r;
    gsl_ieee_double_to_rep (&d, &r);

    gsl_test_int (r.sign, 0, "x = DBL_MAX, sign is +");
    gsl_test_int (r.exponent, 1023, "x = DBL_MAX, exponent is 1023");
    gsl_test_str (r.mantissa, mantissa, "x = DBL_MAX, mantissa");
    gsl_test_int (r.type, GSL_IEEE_TYPE_NORMAL,
		  "x = DBL_MAX, type is NORMAL");
  }

  /* Check for DENORMAL numbers (e.g. DBL_MIN/2^n) */

  {
    double d = 2.2250738585072014e-308;		/* DBL_MIN */
    char mantissa[]
    = "1000000000000000000000000000000000000000000000000000";
    int i;
    gsl_ieee_double_rep r;

    for (i = 0; i < 52; i++)
      {
	double x = d / pow (2.0, 1 + (double) i);
	mantissa[i] = '1';
	gsl_ieee_double_to_rep (&x, &r);

	gsl_test_int (r.sign, 0, "x = DBL_MIN/2^%d, sign is +", i + 1);
	gsl_test_int (r.exponent, -1023,
		      "x = DBL_MIN/2^%d, exponent is -1022", i + 1);
	gsl_test_str (r.mantissa, mantissa,
		      "x = DBL_MIN/2^%d, mantissa", i + 1);
	gsl_test_int (r.type, GSL_IEEE_TYPE_DENORMAL,
		      "x = DBL_MIN/2^%d, type is DENORMAL", i + 1);
	mantissa[i] = '0';
      }
  }

  /* Check for positive INFINITY (e.g. 2*DBL_MAX) */

  {
    double d = 1.7976931348623157e+308;		/* DBL_MAX */
    const char mantissa[]
    = "0000000000000000000000000000000000000000000000000000";
    gsl_ieee_double_rep r;

    double x = 2 * d;
    gsl_ieee_double_to_rep (&x, &r);

    gsl_test_int (r.sign, 0, "x = 2*DBL_MAX, sign is +");
    gsl_test_int (r.exponent, 1024, "x = 2*DBL_MAX, exponent is 1024");
    gsl_test_str (r.mantissa, mantissa, "x = 2*DBL_MAX, mantissa");
    gsl_test_int (r.type, GSL_IEEE_TYPE_INF, "x = 2*DBL_MAX, type is INF");
  }

  /* Check for negative INFINITY (e.g. -2*DBL_MAX) */

  {
    double d = 1.7976931348623157e+308;		/* DBL_MAX */
    const char mantissa[]
    = "0000000000000000000000000000000000000000000000000000";
    gsl_ieee_double_rep r;

    double x = -2 * d;
    gsl_ieee_double_to_rep (&x, &r);

    gsl_test_int (r.sign, 1, "x = -2*DBL_MAX, sign is -");
    gsl_test_int (r.exponent, 1024, "x = -2*DBL_MAX, exponent is 1024");
    gsl_test_str (r.mantissa, mantissa, "x = -2*DBL_MAX, mantissa");
    gsl_test_int (r.type, GSL_IEEE_TYPE_INF,"x = -2*DBL_MAX, type is INF");
  }

  /* Check for NAN (e.g. Inf - Inf) */

  {
    const char mantissa[]
    = "0000000000000000000000000000000000000000000000000000";
    gsl_ieee_double_rep r;
    double x, y, z;

    x = 1.0 / 0.0;
    y = 2.0 / 0.0;
    z = y - x;

    gsl_ieee_double_to_rep (&z, &r);

    /* We don't check the sign and we don't check the mantissa because
       they could be anything for a NaN */

    gsl_test_int (r.exponent, 1024, "x = NaN, exponent is 1024");
    gsl_test_int (r.type, GSL_IEEE_TYPE_NAN, "x = NaN, type is NAN");
  }

  return gsl_test_summary ();
}
