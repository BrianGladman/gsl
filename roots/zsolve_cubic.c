/* zsolve_cubic.c - finds the complex roots of x^3 + a x^2 + b x + c = 0 */

#include <config.h>
#include <math.h>
#include <gsl_math.h>
#include <gsl_complex.h>
#include <gsl_roots.h>

#define SWAP(a,b) do { double tmp = b ; b = a ; a = tmp ; } while(0)

int
gsl_root_complex_solve_cubic (double a, double b, double c, gsl_complex z[])
{
  double Q = (a * a - 3 * b) / 9;
  double R = (2 * a * a * a - 9 * a * b + 27 * c) / 54;
  double Q3 = Q * Q * Q;
  double R2 = R * R;

  if (R == 0 && Q == 0)
    {
      GSL_REAL (z[0]) = -a / 3;
      GSL_IMAG (z[0]) = 0;
      GSL_REAL (z[1]) = -a / 3;
      GSL_IMAG (z[1]) = 0;
      GSL_REAL (z[2]) = -a / 3;
      GSL_IMAG (z[2]) = 0;
      return 3;
    }
  else if (R2 == Q3)
    {

      /* Due to finite precision some double roots may be missed, and
         considered to be a pair of complex roots z = x +/- epsilon i
         close to the real axis. */

      double sqrtQ = sqrt (Q);
      if (R > 0)
	{
	  GSL_REAL (z[0]) = -2 * sqrtQ - a / 3;
	  GSL_IMAG (z[0]) = 0;
	  GSL_REAL (z[1]) = sqrtQ - a / 3;
	  GSL_IMAG (z[1]) = 0;
	  GSL_REAL (z[2]) = sqrtQ - a / 3;
	  GSL_IMAG (z[2]) = 0;
	}
      else
	{
	  GSL_REAL (z[0]) = -sqrtQ - a / 3;
	  GSL_IMAG (z[0]) = 0;
	  GSL_REAL (z[1]) = -sqrtQ - a / 3;
	  GSL_IMAG (z[1]) = 0;
	  GSL_REAL (z[2]) = 2 * sqrtQ - a / 3;
	  GSL_IMAG (z[2]) = 0;
	}
      return 3;
    }
  else if (R2 < Q3)
    {
      double sqrtQ = sqrt (Q);
      double sqrtQ3 = sqrtQ * sqrtQ * sqrtQ;
      double theta = acos (R / sqrtQ3);
      double norm = -2 * sqrtQ;
      double r0 = norm * cos (theta / 3) - a / 3;
      double r1 = norm * cos ((theta + 2.0 * M_PI) / 3) - a / 3;
      double r2 = norm * cos ((theta - 2.0 * M_PI) / 3) - a / 3;

      /* Sort r0, r1, r2 into increasing order */

      if (r0 > r1)
	SWAP (r0, r1);

      if (r1 > r2)
	{
	  SWAP (r1, r2);

	  if (r0 > r1)
	    SWAP (r0, r1);
	}

      GSL_REAL (z[0]) = r0;
      GSL_IMAG (z[0]) = 0;

      GSL_REAL (z[1]) = r1;
      GSL_IMAG (z[1]) = 0;

      GSL_REAL (z[2]) = r2;
      GSL_IMAG (z[2]) = 0;

      return 3;
    }
  else
    {
      double sgnR = (R >= 0 ? 1 : -1);
      double A = -sgnR * pow (fabs (R) + sqrt (R2 - Q3), 1.0 / 3.0);
      double B = Q / A;

      if (A + B < 0)
	{
	  GSL_REAL (z[0]) = A + B - a / 3;
	  GSL_IMAG (z[0]) = 0;

	  GSL_REAL (z[1]) = -0.5 * (A + B) - a / 3;
	  GSL_IMAG (z[1]) = -(sqrt (3.0) / 2.0) * fabs(A - B);

	  GSL_REAL (z[2]) = -0.5 * (A + B) - a / 3;
	  GSL_IMAG (z[2]) = (sqrt (3.0) / 2.0) * fabs(A - B);
	}
      else
	{
	  GSL_REAL (z[0]) = -0.5 * (A + B) - a / 3;
	  GSL_IMAG (z[0]) = -(sqrt (3.0) / 2.0) * fabs(A - B);

	  GSL_REAL (z[1]) = -0.5 * (A + B) - a / 3;
	  GSL_IMAG (z[1]) = (sqrt (3.0) / 2.0) * fabs(A - B);

	  GSL_REAL (z[2]) = A + B - a / 3;
	  GSL_IMAG (z[2]) = 0;
	}

      return 3;
    }
}
