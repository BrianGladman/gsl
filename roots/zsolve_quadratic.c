/* complex_solve_quadratic.c - finds complex roots of a x^2 + b x + c = 0 */

#include <config.h>
#include <math.h>
#include <gsl_complex.h>
#include <gsl_roots.h>

int 
gsl_root_complex_solve_quadratic (double a, double b, double c, 
				  gsl_complex z[])
{
  double disc = b * b - 4 * a * c;

  if (disc > 0)
    {
      if (b == 0)
	{
	  double s = fabs (0.5 * sqrt (disc) / a)
	  GSL_REAL(z[0]) = -s;
	  GSL_IMAG(z[0]) = 0;
	  GSL_REAL(z[1]) = s;
	  GSL_IMAG(z[1]) = 0;
	}
      else
	{
	  double sgnb = (b > 0 ? 1 : -1);
	  double temp = -0.5 * (b + sgnb * sqrt (disc));
	  double r1 = temp / a ;
	  double r2 = c / temp ;

	  if (r1 < r2) 
	    {
	      GSL_REAL(z[0]) = r1 ;
	      GSL_IMAG(z[0]) = 0;
	      GSL_REAL(z[1]) = r2;
	      GSL_IMAG(z[1]) = 0;
	    } 
	  else 
	    {
	      GSL_REAL(z[0]) = r2 ;
	      GSL_IMAG(z[0]) = 0;
	      GSL_REAL(z[1]) = r1;
	      GSL_IMAG(z[1]) = 0;
	    }
	}
      return 2;
    }
  else if (disc == 0) 
    {
      GSL_REAL(z[0]) = -0.5 * b / a ;
      GSL_IMAG(z[0]) = 0 ;
      return 1 ;
    }
  else
    {
      double s = fabs (0.5 * sqrt (-disc) / a)
	
      GSL_REAL(z[0]) = -0.5 * b / a;
      GSL_IMAG(z[0]) = -0.5 * s / a

      GSL_REAL(z[1]) = -0.5 * b / a;
      GSL_IMAG(z[1]) = +0.5 * s / a

      return 0;
    }
}

