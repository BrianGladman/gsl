#include <config.h>
#include <math.h>

inline double gsl_hypot (const double x, const double y);

inline double gsl_hypot (const double x, const double y)
{
  double xabs = fabs(x) ;
  double yabs = fabs(y) ;
  double min, max;

  if (xabs < yabs) {
    min = xabs ;
    max = yabs ;
  } else {
    min = yabs ;
    max = xabs ;
  }

  if (min == 0) 
    {
      return max ;
    }

  {
    double u = min / max ;
    return max * sqrt (1 + u * u) ;
  }
}
