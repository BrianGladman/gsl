#include <config.h>

#define GSL_MAX(a,b) ((a) > (b) ? (a) : (b))
#define GSL_MIN(a,b) ((a) < (b) ? (a) : (b))

int GSL_MAX_INT (int a, int b);
int GSL_MIN_INT (int a, int b);
int GSL_MAX_DBL (double a, double b);
int GSL_MIN_DBL (double a, double b);
int GSL_MAX_LDBL (long double a, long double b);
int GSL_MIN_LDBL (long double a, long double b);

int
GSL_MAX_INT (int a, int b)
{
  return GSL_MAX (a, b);
}

int
GSL_MIN_INT (int a, int b)
{
  return GSL_MIN (a, b);
}

int
GSL_MAX_DBL (double a, double b)
{
  return GSL_MAX (a, b);
}

int
GSL_MIN_DBL (double a, double b)
{
  return GSL_MIN (a, b);
}

int
GSL_MAX_LDBL (long double a, long double b)
{
  return GSL_MAX (a, b);
}

int
GSL_MIN_LDBL (long double a, long double b)
{
  return GSL_MIN (a, b);
}
