#include <stdio.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex_math.h>

struct f1
  {
    char *name;
      gsl_complex (*f) (gsl_complex z);
    double x;
    double y;
    double fx;
    double fy;
  };

#define FN(x) "gsl_complex_" #x, gsl_complex_ ## x
#define ARG(x,y) x, y
#define RES(x,y) x, y

struct f1 list1[] =
{
#include "results.h"
  {"", 0, 0, 0, 0, 0}
};

int
main (void)
{
  size_t i = 0;

  gsl_ieee_env_setup();

  while (list1[i].f)
    {
      struct f1 t = list1[i];
      gsl_complex z = gsl_complex_xy (t.x, t.y);
      gsl_complex fz = (t.f) (z);
      double fx = GSL_REAL (fz), fy = GSL_IMAG (fz);
      gsl_test_rel (fx, t.fx, 10 * GSL_DBL_EPSILON, "%s real part at (%g,%g)", t.name, t.x, t.y);
      gsl_test_rel (fy, t.fy, 10 * GSL_DBL_EPSILON, "%s imag part at (%g,%g)", t.name, t.x, t.y);
      i++;
    }

  return gsl_test_summary ();
}
