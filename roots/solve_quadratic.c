/* solve_quadratic.c - finds the real roots of a x^2 + b x + c = 0 */

#include <config.h>
#include <math.h>

int 
gsl_root_solve_quadratic (double a, double b, double c, double x[])
{
  double disc = b * b - 4 * a * c;

  if (disc > 0)
    {
      if (b == 0)
	{
	  x[0] = -fabs (0.5 * sqrt (disc) / a);
	  x[1] = -x[0];
	}
      else
	{
	  double sgnb = (b > 0 ? 1 : -1);
	  double temp = -0.5 * (b + sgnb * sqrt (disc));
	  double r1 = temp / a ;
	  double r2 = c / temp ;

	  if (r1 < r2) 
	    {
	      x[0] = r1 ;
	      x[1] = r2 ;
	    } 
	  else 
	    {
	      x[0] = r2 ;
	      x[1] = r1 ;
	    }
	}
      return 2;
    }
  else if (disc == 0) 
    {
      x[0] = -0.5 * b / a ;
      return 1 ;
    }
  else
    {
      return 0;
    }
}
