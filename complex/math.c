#include <config.h>
#include <math.h>
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

  double r = gsl_complex_abs (a);
  double theta = gsl_complex_arg (a);

  if (r == 0.0)
    {
      GSL_SET_COMPLEX (&z, 0.0, 0.0);
    }
  else
    {
      double ar = GSL_REAL (a), ai = GSL_IMAG (a);
      double br = GSL_REAL (b), bi = GSL_IMAG (b);

      double rho = pow (r, br) * exp (-bi * theta);
      double beta = theta * br + bi * log (r);

      GSL_SET_COMPLEX (&z, rho * cos (beta), rho * sin (beta));
    }
  
  return z;
}

gsl_complex
gsl_complex_pow_real (gsl_complex a, double b)
{				/* z=a^b */
  gsl_complex z;

  double r = gsl_complex_abs (a);
  double theta = gsl_complex_arg (a);

  if (r == 0)
    {
      GSL_SET_COMPLEX (&z, 0, 0);
    } 
  else 
    {
      double rho = pow (r, b);
      double beta = theta * b;
      GSL_SET_COMPLEX (&z, rho * cos (beta), rho * sin (beta));
    }

  return z;
}

gsl_complex
gsl_complex_log (gsl_complex a)
{				/* z=log(a) */
  double r = gsl_complex_abs (a);
  double theta = gsl_complex_arg (a);

  gsl_complex z;
  GSL_SET_COMPLEX (&z, log(r), theta);
  return z;
}

gsl_complex
gsl_complex_log10 (gsl_complex a)
{				/* r=log10(a) */
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
{				/* r=sin(a) */
  double R = GSL_REAL (a), I = GSL_IMAG (a);
  
  gsl_complex z;
  GSL_SET_COMPLEX (&z, sin (R) * cosh (I), cos (R) * sinh (I));
  return z;
}

gsl_complex
gsl_complex_cos (gsl_complex a)
{				/* r=cos(a) */
  double R = GSL_REAL (a), I = GSL_IMAG (a);
  
  gsl_complex z;
  GSL_SET_COMPLEX (&z, cos (R) * cosh (I), -sin (R) * sinh (I));
  return z;
}

gsl_complex
gsl_complex_sec (gsl_complex a)
{				/* r=sec(a) */
  double R, I, D;
  R = GSL_REAL (a), I = GSL_IMAG (a);
  D = 1 - pow (sin (R), 2.0) + pow (sinh (I), 2.0);
  GSL_SET_COMPLEX (r, cos (R) * cosh (I) / D, sin (R) * sinh (I) / D);
}

gsl_complex
gsl_complex_csc (gsl_complex a)
{				/* r=csc(a) */
  double R, I, D;
  R = GSL_REAL (a), I = GSL_IMAG (a);
  D = pow (sin (R), 2.0) + pow (sinh (I), 2.0);
  GSL_SET_COMPLEX (r, sin (R) * cosh (I) / D, -cos (R) * sinh (I) / D);
}

gsl_complex
gsl_complex_tan (gsl_complex a)
{				/* r=tan(a) */
  double R, I, D;
  R = GSL_REAL (a), I = GSL_IMAG (a);
  D = pow (cos (R), 2.0) + pow (sinh (I), 2.0);
  GSL_SET_COMPLEX (r, 0.5 * sin (2 * R) / D, 0.5 * sinh (2 * I) / D);
}

gsl_complex
gsl_complex_cot (gsl_complex a)
{				/* r=cot(a) */
  double R, I, D;
  R = GSL_REAL (a), I = GSL_IMAG (a);
  D = pow (sin (R), 2.0) + pow (sinh (I), 2.0);
  GSL_SET_COMPLEX (r, 0.5 * sin (2 * R) / D, -0.5 * sinh (2 * I) / D);
}

/**********************************************************************
 * Inverse Complex Trigonometric Functions 
 **********************************************************************/

gsl_complex
gsl_complex_arcsin (gsl_complex a)
{				/* r=arcsin(a) */
  double R, I, D, Phi;
  R = GSL_REAL (a), I = GSL_IMAG (a);
  D = gsl_complex_abs2 (a);
  Phi = sqrt (0.5 * (1 + D + sqrt (4 * pow (I, 2.0) + pow (D - 1, 2.0))));
  gsl_complex_arccosh_d (r, Phi);
  D = GSL_REAL (r);
  gsl_complex_arcsin_d (r, R / Phi);
  GSL_SET_IMAG (r, (I < 0) ? -D : D);
}

gsl_complex
gsl_complex_arcsin_d (double a)
{				/* r=arcsin(a) */
  if (fabs (a) <= 1.0)
    {
      GSL_SET_COMPLEX (r, asin (a), 0.0);
    }
  else
    {
      if (a < 0.0)
	{
	  GSL_SET_COMPLEX (r, -M_PI_2, acosh (-a));
	}
      else
	{
	  GSL_SET_COMPLEX (r, M_PI_2, -acosh (a));
	}
    }
}

gsl_complex
gsl_complex_arccos (gsl_complex a)
{				/* r=arccos(a) */
  double R, I, D, Phi;
  R = GSL_REAL (a), I = GSL_IMAG (a);
  D = gsl_complex_abs2 (a);
  Phi = sqrt (0.5 * (1 + D + sqrt (4 * pow (I, 2.0) + pow (D - 1, 2.0))));
  gsl_complex_arccosh_d (r, Phi);
  D = GSL_REAL (r);
  gsl_complex_arccos_d (r, R / Phi);
  GSL_SET_IMAG (r, (I < 0) ? D : -D);
}

gsl_complex
gsl_complex_arccos_d (double a)
{				/* r=arccos(a) */
  if (fabs (a) <= 1.0)
    {
      GSL_SET_COMPLEX (r, acos (a), 0);
    }
  else
    {
      if (a < 0.0)
	{
	  GSL_SET_COMPLEX (r, M_PI, -acosh (-a));
	}
      else
	{
	  GSL_SET_COMPLEX (r, 0, acosh (a));
	}
    }
}

gsl_complex
gsl_complex_arcsec (gsl_complex a)
{				/* r=arcsec(a) */
  double R, I, D, Phi, P;
  R = GSL_REAL (a), I = GSL_IMAG (a);
  D = pow (R, 2.0) + pow (I, 2.0);
  Phi = sqrt (pow (pow (R, 2.0) - 1, 2.0) + pow (I, 2.0) * (2 * (1 + pow (R, 2.0)) + pow (I, 2.0)));
  Phi = sqrt ((1 + D + Phi) / (2 * D));
  P = R / D;
  GSL_SET_COMPLEX (r, acos (P / Phi), (I >= 0.0) ? acosh (Phi) : -acosh (Phi));
}
gsl_complex
gsl_complex_arcsec_d (double a)
{				/* r=arcsec(a) */
  if (a <= -1.0 || a >= 1.0)
    {
      GSL_SET_COMPLEX (r, asec (a), 0.0);
    }
  else
    {
      if (a >= 0.0)
	{
	  GSL_SET_COMPLEX (r, 0, acosh (1 / a));
	}
      else
	{
	  GSL_SET_COMPLEX (r, M_PI, -acosh (-1 / a));
	}
    }
}

gsl_complex
gsl_complex_arccsc (gsl_complex a)
{				/* r=arccsc(a) */
  double R, I, D, Phi, P;
  R = GSL_REAL (a), I = GSL_IMAG (a);
  D = pow (R, 2.0) + pow (I, 2.0);
  Phi = sqrt (pow (pow (R, 2.0) - 1, 2.0) + pow (I, 2.0) * (2 * (1 + pow (R, 2.0)) + pow (I, 2.0)));
  Phi = sqrt ((1 + D + Phi) / (2 * D));
  P = R / D;
  GSL_SET_COMPLEX (r, asin (P / Phi), (I >= 0.0) ? -acosh (Phi) : acosh (Phi));
}

gsl_complex
gsl_complex_arccsc_d (double a)
{				/* r=arccsc(a) */
  if (a <= -1.0 || a >= 1.0)
    {
      GSL_SET_COMPLEX (r, acsc (a), 0.0);
    }
  else
    {
      if (a >= 0.0)
	{
	  GSL_SET_COMPLEX (r, M_PI_2, -acosh (1 / a));
	}
      else
	{
	  GSL_SET_COMPLEX (r, -M_PI_2, -acosh (-1 / a));
	}
    }
}

gsl_complex
gsl_complex_arctan (gsl_complex a)
{				/* r=arctan(a) */
  double x, R = GSL_REAL (a), I = GSL_IMAG (a);
  GSL_SET_COMPLEX (r, -I, R);	/* *=I */
  gsl_complex_arctanh (r, r);
  x = GSL_REAL (r);
  R = GSL_IMAG (r);
  GSL_SET_COMPLEX (r, R, -x);	/* /=I */
}

gsl_complex
gsl_complex_arccot (gsl_complex a)
{				/* r=arccot(a) */
  if (GSL_REAL (a) == 0.0 && GSL_IMAG (a) == 0.0)
    {
      GSL_SET_COMPLEX (r, M_PI_2, 0);
    }
  else
    {
      gsl_complex_inverse (r, a);
      gsl_complex_arctan (r, r);
    }
}

/**********************************************************************
 * Complex Hyperbolic Functions 
 **********************************************************************/

gsl_complex
gsl_complex_sinh (gsl_complex a)
{				/* r=sinh(a) */
  double R, I;
  R = GSL_REAL (a), I = GSL_IMAG (a);
  GSL_SET_COMPLEX (r, sinh (R) * cos (I), cosh (R) * sin (I));
}

gsl_complex
gsl_complex_cosh (gsl_complex a)
{				/* r=cosh(a) */
  double R, I;
  R = GSL_REAL (a), I = GSL_IMAG (a);
  GSL_SET_COMPLEX (r, cosh (R) * cos (I), sinh (R) * sin (I));
}

gsl_complex
gsl_complex_sech (gsl_complex a)
{				/* r=sech(a) */
  double R, I, D;
  R = GSL_REAL (a), I = GSL_IMAG (a);
  D = pow (cos (I), 2.0) + pow (sinh (R), 2.0);
  GSL_SET_COMPLEX (r, cosh (R) * cos (I) / D, -sinh (R) * sin (I) / D);
}

gsl_complex
gsl_complex_csch (gsl_complex a)
{				/* r=csch(a) */
  double R, I, D;
  R = GSL_REAL (a), I = GSL_IMAG (a);
  D = pow (sinh (R), 2.0) + pow (sin (I), 2.0);
  GSL_SET_COMPLEX (r, sinh (R) * cos (I) / D, -cosh (R) * sin (I) / D);
}

gsl_complex
gsl_complex_tanh (gsl_complex a)
{				/* r=tanh(a) */
  double R, I, D;
  R = GSL_REAL (a), I = GSL_IMAG (a);
  D = pow (cos (I), 2.0) + pow (sinh (R), 2.0);
  GSL_SET_COMPLEX (r, sinh (R) * cosh (R) / D, 0.5 * sin (2 * I) / D);
}

gsl_complex
gsl_complex_coth (gsl_complex a)
{				/* r=coth(a) */
  double R, I, D;
  R = GSL_REAL (a), I = GSL_IMAG (a);
  D = pow (sin (I), 2.0) + pow (sinh (R), 2.0);
  GSL_SET_COMPLEX (r, sinh (R) * cosh (R) / D, -0.5 * sin (2 * I) / D);
}

/**********************************************************************
 * Inverse Complex Hyperbolic Functions 
 **********************************************************************/

gsl_complex
gsl_complex_arcsinh (gsl_complex a)
{				/* r=arcsinh(a) */
  double R, I, T;
  R = GSL_REAL (a), I = GSL_IMAG (a);
  if (R == 0)
    {
      if (I == 0)
	{
	  GSL_SET_COMPLEX (r, 0.0, 0.0);
	}
      else
	{
	  T = 1 - pow (I, 2.0);
	  if (T >= 0)
	    {
	      GSL_SET_COMPLEX (r, 0, atan2 (I, sqrt (T)));
	    }
	  else
	    {
	      T = sqrt (-T);
	      if (I >= 0)
		T += I;
	      else
		T -= I;
	      if (I >= 0)
		{
		  GSL_SET_COMPLEX (r, log (T), M_PI_2);
		}
	      else
		{
		  GSL_SET_COMPLEX (r, -log (T), -M_PI_2);
		}
	    }
	}
    }
  else
    {
      if (I == 0)
	{
	  T = sqrt (1 + pow (R, 2.0));
	  if (fabs (R) < 0.5)
	    {
	      GSL_SET_COMPLEX (r, atanh (R / T), 0);
	    }
	  else
	    {
	      if (R >= 0.5)
		GSL_SET_COMPLEX (r, log (T + R), 0);
	      else
		GSL_SET_COMPLEX (r, -log (T - R), 0);
	    }
	}
      else
	{
	  gsl_complex z;
	  gsl_complex_mul (&z, a, a);	/* z=a^2 */
	  gsl_complex_add_d (&z, &z, 1.0);
	  gsl_complex_sqrt (&z, &z);	/* z=sqrt(a^2+1) */
	  gsl_complex_add_d (&z, &z, 1.0);
	  gsl_complex_div (&z, a, &z);	/* z=a/(1+sqrt(a^2+1) */
	  gsl_complex_arctanh (r, &z);
	  gsl_complex_mul_d (r, r, 2.0);
	}
    }
}

gsl_complex
gsl_complex_arccosh (gsl_complex a)
{				/* r=arccosh(a) */
  double R, I;
  gsl_complex z;
  R = GSL_REAL (a), I = GSL_IMAG (a);
  if (I == 0)
    {
      if (R == 0)
	{
	  GSL_SET_COMPLEX (r, 0, M_PI_2);
	  return;
	}
      if (R == 1)
	{
	  GSL_SET_COMPLEX (r, 0.0, 0.0);
	  return;
	}
      if (R == -1)
	{
	  GSL_SET_COMPLEX (r, 0, M_PI);
	  return;
	}
      if (R < -1)
	{
	  GSL_SET_COMPLEX (r, log (sqrt (pow (R, 2.0) - 1) - R), M_PI);
	  return;
	}
    }
  gsl_complex_add_d (&z, a, 1.0);	/* z=1+a */
  gsl_complex_sub_d (r, a, 1.0);	/* r=a-1 */
  gsl_complex_mul_d (r, r, 0.5);
  gsl_complex_sqrt (r, r);	/* r=sqrt(0.5*(1+a)); */
  gsl_complex_mul_d (&z, &z, 0.5);
  gsl_complex_sqrt (&z, &z);
  gsl_complex_add_d (&z, &z, 1.0);	/* z=1+sqrt(0.5*(1+a)) */
  gsl_complex_div (r, r, &z);	/* r=sqrt(0.5*(1+a))/(1+sqrt(0.5*(1+a))) */
  gsl_complex_arctanh (r, r);
  gsl_complex_mul_d (r, r, 4.0);
}

gsl_complex
gsl_complex_arccosh_d (double a)
{				/* r=arccosh(a) */
  if (a >= 1)
    {
      GSL_SET_COMPLEX (r, acosh (a), 0);
    }
  else
    {
      if (a >= -1.0)
	{
	  GSL_SET_COMPLEX (r, 0, acos (a));
	}
      else
	{
	  GSL_SET_COMPLEX (r, acosh (-a), M_PI);
	}
    }
}

gsl_complex
gsl_complex_arcsech (gsl_complex a)
{				/* r=arcsech(a); */
  gsl_complex_inverse (r, a);
  gsl_complex_arccosh (r, r);
}

gsl_complex
gsl_complex_arccsch (gsl_complex a)
{				/* r=arccsch(a) */
  gsl_complex_inverse (r, a);
  gsl_complex_arcsinh (r, r);
}

gsl_complex
gsl_complex_arctanh (gsl_complex a)
{				/* r=arctanh(a) */
  double R, I, A, B, X, Y, T1, T2, D;
  R = GSL_REAL (a), I = GSL_IMAG (a);
  if (R == 0.0)
    {
      GSL_SET_COMPLEX (r, 0, atan (I));
    }
  else
    {
      if (I == 0.0)
	{
	  if (fabs (R) < 0.5)
	    {
	      GSL_SET_COMPLEX (r, atanh (R), 0.0);
	    }
	  else
	    {
	      A = 1 - R;
	      B = (1 + R) / A;
	      if (B >= 0)
		{
		  GSL_SET_IMAG (r, 0.0);
		}
	      else
		{
		  B = -B;
		  GSL_SET_IMAG (r, (A < 0) ? -M_PI_2 : M_PI_2);
		}
	      GSL_SET_REAL (r, 0.5 * log (B));
	    }
	}
      else
	{
	  D = gsl_complex_abs2 (a);
	  A = 1 + R, B = 1 - R;
	  T1 = fabs (4 * R), T2 = 1 + D;
	  if (T1 < T2)
	    {
	      GSL_SET_REAL (r, 0.5 * atanh (2 * R / T2));
	    }
	  else
	    {
	      GSL_SET_REAL (r, 0.25 * log ((pow (A, 2.0) + pow (I, 2.0)) / (pow (B, 2.0) + pow (I, 2.0))));
	    }
	  X = A * B - pow (I, 2.0);
	  Y = 2 * I;
	  GSL_SET_IMAG (r, 0.5 * atan2 (Y, X));
	}
    }
}

gsl_complex
gsl_complex_arctanh_d (double a)
{				/* r=arctanh(a) */
  if (a > -1.0 && a < 1.0)
    {
      GSL_SET_COMPLEX (r, atanh (a), 0);
    }
  else
    {
      GSL_SET_COMPLEX (r, atanh (1 / a), (a < 0) ? M_PI_2 : -M_PI_2);
    }
}

gsl_complex
gsl_complex_arccoth (gsl_complex a)
{				/* r=arccoth(a) */
  gsl_complex_inverse (r, a);
  gsl_complex_arctanh (r, r);
}
