/* solve_cubic.c - finds the real roots of x^3 + a x^2 + b x + c = 0 */

#include <config.h>
#include <math.h>
#include <gsl_math.h>
#include <gsl_poly.h>

#define SWAP(a,b) do { double tmp = b ; b = a ; a = tmp ; } while(0)

int 
gsl_poly_solve_cubic (double a, double b, double c, double x[])
{
  double Q = (a * a - 3 * b) / 9;
  double R = (2 * a * a * a - 9 * a * b + 27 * c) / 54;
  double Q3 = Q * Q * Q;
  double R2 = R * R;

  if (R == 0 && Q == 0)
    {
      x[0] = - a / 3 ;
      x[1] = - a / 3 ;
      x[2] = - a / 3 ;
      return 3 ;
    }
  else if (R2 == Q3)
    {

      /* Due to finite precision some double roots may be missed, and
	 considered to be a pair of complex roots z = x +/- epsilon i
	 close to the real axis. */

      double sqrtQ = sqrt (Q);
      if (R > 0)
	{
	  x[0] = -2 * sqrtQ  - a / 3;
	  x[1] = sqrtQ - a / 3;
	  x[2] = sqrtQ - a / 3;
	}
      else
	{
	  x[0] = - sqrtQ  - a / 3;
	  x[1] = - sqrtQ - a / 3;
	  x[2] = 2 * sqrtQ - a / 3;
	}
      return 3 ;
    }
  else if (R2 < Q3)
    {
      double sqrtQ = sqrt (Q);
      double sqrtQ3 = sqrtQ * sqrtQ * sqrtQ;
      double theta = acos (R / sqrtQ3);
      double norm = -2 * sqrtQ;
      x[0] = norm * cos (theta / 3) - a / 3;
      x[1] = norm * cos ((theta + 2.0 * M_PI) / 3) - a / 3;
      x[2] = norm * cos ((theta - 2.0 * M_PI) / 3) - a / 3;
      
      /* Sort x[0], x[1], x[2] into increasing order */

      if (x[0] > x[1])
	SWAP(x[0], x[1]) ;
      
      if (x[1] > x[2])
	{
	  SWAP(x[1], x[2]) ;
	  
	  if (x[0] > x[1])
	    SWAP(x[0], x[1]) ;
	}
      
      return 3;
    }
  else
    {
      double sgnR = (R >= 0 ? 1 : -1);
      double A = -sgnR * pow (fabs (R) + sqrt (R2 - Q3), 1.0/3.0);
      double B = Q / A ;
      x[0] = A + B - a / 3;
      return 1;
    }
}
