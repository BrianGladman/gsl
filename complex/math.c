/* Basic complex arithmetic functions 
 *
 * Original version by Jorma Olavi T{htinen <jotahtin@cc.hut.fi> 
 *
 * Modified for GSL by Brian Gough, 3/2000 
 */

#include <config.h>
#include <math.h>
#include <gsl_math.h>
#include <gsl_complex.h>
#include <gsl_complex_math.h>

/**********************************************************************
 * Complex numbers 
 **********************************************************************/

gsl_complex
gsl_complex_xy (double x, double y)
{				/* return z = x + i y */
  gsl_complex z;
  GSL_SET_COMPLEX (&z, x, y);
  return z;
}

gsl_complex
gsl_complex_polar (double r, double theta)
{				/* return z = r exp(i theta) */
  gsl_complex z;
  GSL_SET_COMPLEX (&z, r * cos (theta), r * sin (theta));
  return z;
}

/**********************************************************************
 * Properties of complex numbers 
 **********************************************************************/

double
gsl_complex_arg (gsl_complex z)
{				/* return arg(z),  -pi < arg(z) <= +pi */
  return atan2 (GSL_IMAG (z), GSL_REAL (z));
}

double
gsl_complex_abs (gsl_complex z)
{				/* return |z| */
  return gsl_hypot (GSL_REAL (z), GSL_IMAG (z));
}

double
gsl_complex_abs2 (gsl_complex z)
{				/* return |z|^2 */
  double x = GSL_REAL (z);
  double y = GSL_IMAG (z);

  return (x * x + y * y);
}

/* The behavior of these functions can be improved by handling
   underflows more gracefully, as shown here: */

double
gsl_complex_logabs (gsl_complex z)
{				/* return log|z| */
  double xabs = fabs(GSL_REAL (z));
  double yabs = fabs(GSL_IMAG (z));
  double max, r;

  if (xabs >= yabs)
    {
      max = xabs;
      r = yabs/xabs;
    }
  else
    {
      max = yabs;
      r = xabs/yabs;
    }

  return log (max) + 0.5 * log1p(r*r) ;   /* FIXME: non-ansi log1p */
}


/***********************************************************************
 * Complex arithmetic operators 
 ***********************************************************************/

gsl_complex
gsl_complex_add (gsl_complex a, gsl_complex b)
{				/* z=a+b */
  double ar = GSL_REAL (a), ai = GSL_IMAG (a);
  double br = GSL_REAL (b), bi = GSL_IMAG (b);

  gsl_complex z;
  GSL_SET_COMPLEX (&z, ar + br, ai + bi);
  return z;
}

gsl_complex
gsl_complex_add_real (gsl_complex a, double b)
{				/* z=a+b */
  gsl_complex z;
  GSL_SET_COMPLEX (&z, GSL_REAL (a) + b, GSL_IMAG (a));
  return z;
}

gsl_complex
gsl_complex_sub (gsl_complex a, gsl_complex b)
{				/* z=a-b */
  double ar = GSL_REAL (a), ai = GSL_IMAG (a);
  double br = GSL_REAL (b), bi = GSL_IMAG (b);

  gsl_complex z;
  GSL_SET_COMPLEX (&z, ar - br, ai - bi);
  return z;
}

gsl_complex
gsl_complex_sub_real (gsl_complex a, double b)
{				/* z=a-b */
  gsl_complex z;
  GSL_SET_COMPLEX (&z, GSL_REAL (a) - b, GSL_IMAG (a));
  return z;
}

gsl_complex
gsl_complex_mul (gsl_complex a, gsl_complex b)
{				/* z=a*b */
  double ar = GSL_REAL (a), ai = GSL_IMAG (a);
  double br = GSL_REAL (b), bi = GSL_IMAG (b);

  gsl_complex z;
  GSL_SET_COMPLEX (&z, ar * br - ai * bi, ar * bi + ai * br);
  return z;
}

gsl_complex
gsl_complex_mul_real (gsl_complex a, double b)
{				/* z=a*b */
  gsl_complex z;
  GSL_SET_COMPLEX (&z, GSL_REAL (a) * b, GSL_IMAG (a) * b);
  return z;
}

gsl_complex
gsl_complex_div (gsl_complex a, gsl_complex b)
{				/* z=a/b */
  double ar = GSL_REAL (a), ai = GSL_IMAG (a);
  double br = GSL_REAL (b), bi = GSL_IMAG (b);

  double s = 1.0 / gsl_complex_abs (b);

  double sbr = s * br;
  double sbi = s * bi;

  double zr = (ar * sbr + ai * sbi) * s;
  double zi = (ai * sbr - ar * sbi) * s;

  gsl_complex z;
  GSL_SET_COMPLEX (&z, zr, zi);
  return z;
}

gsl_complex
gsl_complex_conjugate (gsl_complex a)
{				/* z=conj(a) */
  gsl_complex z;
  GSL_SET_COMPLEX (&z, GSL_REAL (a), -GSL_IMAG (a));
  return z;
}

gsl_complex
gsl_complex_negative (gsl_complex a)
{				/* z=-a */
  gsl_complex z;
  GSL_SET_COMPLEX (&z, -GSL_REAL (a), -GSL_IMAG (a));
  return z;
}

gsl_complex
gsl_complex_inverse (gsl_complex a)
{				/* z=1/a */
  double s2 = 1.0 / gsl_complex_abs2 (a);

  gsl_complex z;
  GSL_SET_COMPLEX (&z, GSL_REAL(a) * s2, -GSL_IMAG(a) * s2);
  return z;
}

/* Elementary complex functions 
 */

gsl_complex
gsl_complex_sqrt (gsl_complex a)
{				/* z=sqrt(a) */
  gsl_complex z;

  if (GSL_REAL (a) == 0.0 && GSL_IMAG (a) == 0.0)
    {
      GSL_SET_COMPLEX (&z, 0, 0);
    }
  else
    {
      double x = fabs (GSL_REAL (a));
      double y = fabs (GSL_IMAG (a));
      double w;

      if (x >= y)
	{
	  double t = y / x;
	  w = sqrt (x) * sqrt (0.5 * (1.0 + sqrt (1.0 + t * t)));
	}
      else
	{
	  double t = x / y;
	  w = sqrt (y) * sqrt (0.5 * (t + sqrt (1.0 + t * t)));
	}

      if (GSL_REAL (a) >= 0.0)
	{
	  double ai = GSL_IMAG (a);
	  GSL_SET_COMPLEX (&z, w, ai / (2.0 * w));
	}
      else
	{
	  double ai = GSL_IMAG (a);
	  double vi = (ai >= 0) ? w : -w;
	  GSL_SET_COMPLEX (&z, ai / (2.0 * vi), vi);
	}
    }

  return z;
}

gsl_complex
gsl_complex_sqrt_real (double x)
{				/* z=sqrt(x) */
  gsl_complex z;

  if (x >= 0)
    {
      GSL_SET_COMPLEX (&z, sqrt (x), 0.0);
    }
  else
    {
      GSL_SET_COMPLEX (&z, 0.0, sqrt (-x));
    }

  return z;
}

gsl_complex
gsl_complex_exp (gsl_complex a)
{				/* z=exp(a) */
  double rho = exp (GSL_REAL (a));
  double theta = GSL_IMAG (a);

  gsl_complex z;
  GSL_SET_COMPLEX (&z, rho * cos (theta), rho * sin (theta));
  return z;
}

gsl_complex
gsl_complex_pow (gsl_complex a, gsl_complex b)
{				/* z=a^b */
  gsl_complex z;

  if (GSL_REAL(a) == 0 && GSL_IMAG(a) == 0.0)
    {
      GSL_SET_COMPLEX (&z, 0.0, 0.0);
    }
  else
    {
      double logr = gsl_complex_logabs (a);
      double theta = gsl_complex_arg (a);

      double br = GSL_REAL (b), bi = GSL_IMAG (b);

      double rho = exp (logr * br - bi * theta);
      double beta = theta * br + bi * logr;

      GSL_SET_COMPLEX (&z, rho * cos (beta), rho * sin (beta));
    }
  
  return z;
}

gsl_complex
gsl_complex_pow_real (gsl_complex a, double b)
{				/* z=a^b */
  gsl_complex z;

  if (GSL_REAL(a) == 0 && GSL_IMAG(a) == 0)
    {
      GSL_SET_COMPLEX (&z, 0, 0);
    } 
  else 
    {
      double logr = gsl_complex_logabs (a);
      double theta = gsl_complex_arg (a);
      double rho = exp (logr *  b);
      double beta = theta * b;
      GSL_SET_COMPLEX (&z, rho * cos (beta), rho * sin (beta));
    }

  return z;
}

gsl_complex
gsl_complex_log (gsl_complex a)
{				/* z=log(a) */
  double logr = gsl_complex_logabs (a);
  double theta = gsl_complex_arg (a);

  gsl_complex z;
  GSL_SET_COMPLEX (&z, logr, theta);
  return z;
}

gsl_complex
gsl_complex_log10 (gsl_complex a)
{				/* z = log10(a) */
  return gsl_complex_mul_real(gsl_complex_log (a), 1/log(10.));
}

gsl_complex
gsl_complex_log_b (gsl_complex a, gsl_complex b)
{
  return gsl_complex_div (gsl_complex_log(a), gsl_complex_log(b));
}

/***********************************************************************
 * Complex trigonometric functions
 ***********************************************************************/

gsl_complex
gsl_complex_sin (gsl_complex a)
{				/* z = sin(a) */
  double R = GSL_REAL (a), I = GSL_IMAG (a);
  
  gsl_complex z;
  GSL_SET_COMPLEX (&z, sin (R) * cosh (I), cos (R) * sinh (I));
  return z;
}

gsl_complex
gsl_complex_cos (gsl_complex a)
{				/* z = cos(a) */
  double R = GSL_REAL (a), I = GSL_IMAG (a);
  
  gsl_complex z;
  GSL_SET_COMPLEX (&z, cos (R) * cosh (I), -sin (R) * sinh (I));
  return z;
}

gsl_complex
gsl_complex_tan (gsl_complex a)
{				/* z = tan(a) */
  double R = GSL_REAL (a), I = GSL_IMAG (a);
  double D = pow (cos (R), 2.0) + pow (sinh (I), 2.0);

  gsl_complex z;
  GSL_SET_COMPLEX (&z, 0.5 * sin (2 * R) / D, 0.5 * sinh (2 * I) / D);
  return z;
}

gsl_complex
gsl_complex_sec (gsl_complex a)
{				/* z = sec(a) */
  double R = GSL_REAL (a), I = GSL_IMAG (a);

  /* Changed the following expression to be purely positive, was
     previously 1 - sin^2 + sinh^2. BJG 22/MAR/00 */

  double D = pow (cos (R), 2.0) + pow (sinh (I), 2.0);  

  gsl_complex z;
  GSL_SET_COMPLEX (&z, cos (R) * cosh (I) / D, sin (R) * sinh (I) / D);
  return z;
}

gsl_complex
gsl_complex_csc (gsl_complex a)
{				/* z = csc(a) */
  double R = GSL_REAL (a), I = GSL_IMAG (a);
  double D = pow (sin (R), 2.0) + pow (sinh (I), 2.0);
  
  gsl_complex z;
  GSL_SET_COMPLEX (&z, sin (R) * cosh (I) / D, -cos (R) * sinh (I) / D);
  return z;
}


gsl_complex
gsl_complex_cot (gsl_complex a)
{				/* z = cot(a) */
  double R = GSL_REAL (a), I = GSL_IMAG (a);
  double D = pow (sin (R), 2.0) + pow (sinh (I), 2.0);
  
  gsl_complex z ;
  GSL_SET_COMPLEX (&z, 0.5 * sin (2 * R) / D, -0.5 * sinh (2 * I) / D);
  return z;
}

/**********************************************************************
 * Inverse Complex Trigonometric Functions 
 **********************************************************************/

gsl_complex
gsl_complex_arcsin (gsl_complex a)
{				/* z = arcsin(a) */
  double R = GSL_REAL (a), I = GSL_IMAG (a);
  double D = gsl_complex_abs2 (a);
  double Phi = sqrt (0.5 * (1 + D + sqrt (4 * pow (I, 2.0) + pow (D - 1, 2.0))));
  gsl_complex z = gsl_complex_arccosh_real (Phi);
  D = GSL_REAL (z);
  z = gsl_complex_arcsin_real (R / Phi);
  GSL_SET_IMAG (&z, (I < 0) ? -D : D);
  return z;
}

gsl_complex
gsl_complex_arcsin_real (double a)
{				/* z = arcsin(a) */
  gsl_complex z;

  if (fabs (a) <= 1.0)
    {
      GSL_SET_COMPLEX (&z, asin (a), 0.0);
    }
  else
    {
      if (a < 0.0)
	{
	  GSL_SET_COMPLEX (&z, -M_PI_2, acosh (-a));
	}
      else
	{
	  GSL_SET_COMPLEX (&z, M_PI_2, -acosh (a));
	}
    }

  return z;
}

gsl_complex
gsl_complex_arccos (gsl_complex a)
{				/* z = arccos(a) */
  double R = GSL_REAL (a), I = GSL_IMAG (a);
  double D = gsl_complex_abs2 (a);
  double Phi = sqrt (0.5 * (1 + D + sqrt (4 * pow (I, 2.0) + pow (D - 1, 2.0))));
  gsl_complex z = gsl_complex_arccosh_real (Phi);
  D = GSL_REAL (z);
  z = gsl_complex_arccos_real (R / Phi);
  GSL_SET_IMAG (&z, (I < 0) ? D : -D);
  return z;
}

gsl_complex
gsl_complex_arccos_real (double a)
{				/* z = arccos(a) */
  gsl_complex z;

  if (fabs (a) <= 1.0)
    {
      GSL_SET_COMPLEX (&z, acos (a), 0);
    }
  else
    {
      if (a < 0.0)
	{
	  GSL_SET_COMPLEX (&z, M_PI, -acosh (-a));
	}
      else
	{
	  GSL_SET_COMPLEX (&z, 0, acosh (a));
	}
    }

  return z;
}

gsl_complex
gsl_complex_arctan (gsl_complex a)
{				/* z = arctan(a) */
  double x, R = GSL_REAL (a), I = GSL_IMAG (a);
  gsl_complex z;
  GSL_SET_COMPLEX (&z, -I, R);	/* *=I */
  z = gsl_complex_arctanh (z);
  x = GSL_REAL (z);
  R = GSL_IMAG (z);
  GSL_SET_COMPLEX (&z, R, -x);	/* /=I */
  return z;
}

gsl_complex
gsl_complex_arcsec (gsl_complex a)
{				/* z = arcsec(a) */
  double R = GSL_REAL (a), I = GSL_IMAG (a);
  double D = pow (R, 2.0) + pow (I, 2.0);
  double t = sqrt (pow (pow (R, 2.0) - 1, 2.0) 
                     + pow (I, 2.0) * (2 * (1 + pow (R, 2.0)) + pow (I, 2.0)));
  double Phi = sqrt ((1 + D + t) / (2 * D));
  double P = R / D;
  gsl_complex z;
  GSL_SET_COMPLEX (&z, acos (P / Phi), (I >= 0.0) ? acosh (Phi) : -acosh (Phi));
  return z;
}
gsl_complex
gsl_complex_arcsec_real (double a)
{				/* z = arcsec(a) */
  gsl_complex z;

  if (a <= -1.0 || a >= 1.0)
    {
      GSL_SET_COMPLEX (&z, acos (1/a), 0.0);
    }
  else
    {
      if (a >= 0.0)
	{
	  GSL_SET_COMPLEX (&z, 0, acosh (1 / a));
	}
      else
	{
	  GSL_SET_COMPLEX (&z, M_PI, -acosh (-1 / a));
	}
    }

  return z;
}

gsl_complex
gsl_complex_arccsc (gsl_complex a)
{				/* z = arccsc(a) */
  double R = GSL_REAL (a), I = GSL_IMAG (a);
  double D = pow (R, 2.0) + pow (I, 2.0);
  double t = sqrt (pow (pow (R, 2.0) - 1, 2.0) + pow (I, 2.0) * (2 * (1 + pow (R, 2.0)) + pow (I, 2.0)));
  double Phi = sqrt ((1 + D + t) / (2 * D));
  double P = R / D;
  gsl_complex z;
  GSL_SET_COMPLEX (&z, asin (P / Phi), (I >= 0.0) ? -acosh (Phi) : acosh (Phi));
  return z;
}

gsl_complex
gsl_complex_arccsc_real (double a)
{				/* z = arccsc(a) */
  gsl_complex z;

  if (a <= -1.0 || a >= 1.0)
    {
      GSL_SET_COMPLEX (&z, asin (1/a), 0.0);
    }
  else
    {
      if (a >= 0.0)
	{
	  GSL_SET_COMPLEX (&z, M_PI_2, -acosh (1 / a));
	}
      else
	{
	  GSL_SET_COMPLEX (&z, -M_PI_2, -acosh (-1 / a));
	}
    }

  return z;
}


gsl_complex
gsl_complex_arccot (gsl_complex a)
{				/* z = arccot(a) */
  gsl_complex z;

  if (GSL_REAL (a) == 0.0 && GSL_IMAG (a) == 0.0)
    {
      GSL_SET_COMPLEX (&z, M_PI_2, 0);
    }
  else
    {
      z = gsl_complex_inverse (a);
      z = gsl_complex_arctan (z);
    }

  return z;
}

/**********************************************************************
 * Complex Hyperbolic Functions 
 **********************************************************************/

gsl_complex
gsl_complex_sinh (gsl_complex a)
{				/* z = sinh(a) */
  double R = GSL_REAL (a), I = GSL_IMAG (a);

  gsl_complex z;
  GSL_SET_COMPLEX (&z, sinh (R) * cos (I), cosh (R) * sin (I));
  return z;
}

gsl_complex
gsl_complex_cosh (gsl_complex a)
{				/* z = cosh(a) */
  double R = GSL_REAL (a), I = GSL_IMAG (a);

  gsl_complex z;
  GSL_SET_COMPLEX (&z, cosh (R) * cos (I), sinh (R) * sin (I));
  return z;
}

gsl_complex
gsl_complex_sech (gsl_complex a)
{				/* z = sech(a) */
  double R = GSL_REAL (a), I = GSL_IMAG (a);
  double D = pow (cos (I), 2.0) + pow (sinh (R), 2.0);

  gsl_complex z;
  GSL_SET_COMPLEX (&z, cosh (R) * cos (I) / D, -sinh (R) * sin (I) / D);
  return z;
}

gsl_complex
gsl_complex_csch (gsl_complex a)
{				/* z = csch(a) */
  double R = GSL_REAL (a), I = GSL_IMAG (a);
  double D = pow (sinh (R), 2.0) + pow (sin (I), 2.0);

  gsl_complex z;
  GSL_SET_COMPLEX (&z, sinh (R) * cos (I) / D, -cosh (R) * sin (I) / D);
  return z;
}

gsl_complex
gsl_complex_tanh (gsl_complex a)
{				/* z = tanh(a) */
  double R = GSL_REAL (a), I = GSL_IMAG (a);
  double D = pow (cos (I), 2.0) + pow (sinh (R), 2.0);

  gsl_complex z;
  GSL_SET_COMPLEX (&z, sinh (R) * cosh (R) / D, 0.5 * sin (2 * I) / D);
  return z;
}

gsl_complex
gsl_complex_coth (gsl_complex a)
{				/* z = coth(a) */
  double R = GSL_REAL (a), I = GSL_IMAG (a);
  double D = pow (sin (I), 2.0) + pow (sinh (R), 2.0);

  gsl_complex z;
  GSL_SET_COMPLEX (&z, sinh (R) * cosh (R) / D, -0.5 * sin (2 * I) / D);
  return z;
}

/**********************************************************************
 * Inverse Complex Hyperbolic Functions 
 **********************************************************************/

gsl_complex
gsl_complex_arcsinh (gsl_complex a)
{				/* z = arcsinh(a) */
  double R = GSL_REAL (a), I = GSL_IMAG (a);

  gsl_complex z;
  GSL_SET_COMPLEX(&z, (R-I)*(R+I) + 1.0, 2 * R * I);
  z = gsl_complex_sqrt(z) ;
  z = gsl_complex_add (z, a);
  return gsl_complex_log(z);
}

gsl_complex
gsl_complex_arccosh (gsl_complex a)
{				/* z = arccosh(a) */
  double R = GSL_REAL (a), I = GSL_IMAG (a);

  gsl_complex z;
  GSL_SET_COMPLEX(&z, (R-I)*(R+I) - 1.0, 2 * R * I);
  z = gsl_complex_sqrt(z) ;
  z = gsl_complex_add (z, a);
  return gsl_complex_log(z);
}

gsl_complex
gsl_complex_arccosh_real (double a)
{				/* z = arccosh(a) */
  gsl_complex z;

  if (a >= 1)
    {
      GSL_SET_COMPLEX (&z, acosh (a), 0);
    }
  else
    {
      if (a >= -1.0)
	{
	  GSL_SET_COMPLEX (&z, 0, acos (a));
	}
      else
	{
	  GSL_SET_COMPLEX (&z, acosh (-a), M_PI);
	}
    }
  
  return z;
}

gsl_complex
gsl_complex_arcsech (gsl_complex a)
{				/* z = arcsech(a); */
  gsl_complex t = gsl_complex_inverse (a);
  return gsl_complex_arccosh (t);
}

gsl_complex
gsl_complex_arccsch (gsl_complex a)
{				/* z = arccsch(a) */
  gsl_complex t = gsl_complex_inverse (a);
  return gsl_complex_arcsinh (t);
}

gsl_complex
gsl_complex_arctanh (gsl_complex a)
{				/* z = arctanh(a) */
  if (GSL_IMAG(a) == 0.0)
    {
      return gsl_complex_arctanh_real (GSL_REAL(a));
    }
  else 
    {
      gsl_complex p = gsl_complex_add_real(a, 1.0);
      gsl_complex m = gsl_complex_add_real(gsl_complex_negative(a), 1.0);
      gsl_complex r = gsl_complex_div (p, m);
      return gsl_complex_mul_real(gsl_complex_log(r), 0.5);
    }
}

gsl_complex
gsl_complex_arctanh_real (double a)
{				/* z = arctanh(a) */
  gsl_complex z;

  if (a > -1.0 && a < 1.0)
    {
      GSL_SET_COMPLEX (&z, atanh (a), 0);
    }
  else
    {
      GSL_SET_COMPLEX (&z, atanh (1 / a), (a < 0) ? M_PI_2 : -M_PI_2);
    }

  return z;
}

gsl_complex
gsl_complex_arccoth (gsl_complex a)
{				/* z = arccoth(a) */
  gsl_complex t = gsl_complex_inverse (a);
  return gsl_complex_arctanh (t);
}

