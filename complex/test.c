#include <stdio.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

struct f
{
  char *name;
  double (*f) (gsl_complex z);
  double x;
  double y;
  double fx;
  double fy;
};

struct fz
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

struct f list[] =
{
#include "results1.h"
  {"", 0, 0, 0, 0, 0}
};


struct fz listz[] =
{
#include "results.h"
  {"", 0, 0, 0, 0, 0}
};

int
main (void)
{
  size_t i = 0;

  gsl_ieee_env_setup();


  for (i = 0 ; i < 10; i++) 
    {
      double r = (i - 5.0) * 0.3 ;
      double t = 2.0 * M_PI * i / 5 ;
      double x = r * cos(t), y = r * sin(t) ;
      gsl_complex z = gsl_complex_polar (r, t) ;
      gsl_test_rel (GSL_REAL(z), x, 10 * GSL_DBL_EPSILON, "gsl_complex_polar real part at (r=%g,t=%g)", r, t);
      
      gsl_test_rel (GSL_IMAG(z), y, 10 * GSL_DBL_EPSILON, "gsl_complex_polar imag part at (r=%g,t=%g)", r, t);
    }
    
    i = 0;

  while (list[i].f)
    {
      struct f t = list[i];
      gsl_complex z = gsl_complex_xy (t.x, t.y);
      double f = (t.f) (z);
      gsl_test_rel (f, t.fx, 10 * GSL_DBL_EPSILON, "%s at (%g,%g)", t.name, t.x, t.y);
      i++;
    }

  i = 0;

  while (listz[i].f)
    {
      struct fz t = listz[i];
      gsl_complex z = gsl_complex_xy (t.x, t.y);
      gsl_complex fz = (t.f) (z);
      double fx = GSL_REAL (fz), fy = GSL_IMAG (fz);
      gsl_test_rel (fx, t.fx, 10 * GSL_DBL_EPSILON, "%s real part at (%g,%g)", t.name, t.x, t.y);
      gsl_test_rel (fy, t.fy, 10 * GSL_DBL_EPSILON, "%s imag part at (%g,%g)", t.name, t.x, t.y);
      i++;
    }

  return gsl_test_summary ();
}
