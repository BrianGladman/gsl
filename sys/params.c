#include <stdio.h>
#include <math.h>
#include <float.h>

#define rt3(x) pow((x), 1.0 / 3.0)
#define rt4(x) pow((x), 1.0 / 4.0)
#define rt5(x) pow((x), 1.0 / 5.0)
#define rt6(x) pow((x), 1.0 / 6.0)

int
main (void)
{
  printf ("#define GSL_DBL_EPSILON       % .16e\n", DBL_EPSILON);
  printf ("#define GSL_SQRT_DBL_EPSILON  % .16e\n", sqrt (DBL_EPSILON));
  printf ("#define GSL_ROOT3_DBL_EPSILON % .16e\n", rt3 (DBL_EPSILON));
  printf ("#define GSL_ROOT4_DBL_EPSILON % .16e\n", rt4 (DBL_EPSILON));
  printf ("#define GSL_ROOT5_DBL_EPSILON % .16e\n", rt5 (DBL_EPSILON));
  printf ("#define GSL_ROOT6_DBL_EPSILON % .16e\n", rt6 (DBL_EPSILON));
  printf ("#define GSL_LOG_DBL_EPSILON   % .16e\n", log (DBL_EPSILON));
  printf ("\n");

  printf ("#define GSL_DBL_MIN       % .16e\n", DBL_MIN);
  printf ("#define GSL_SQRT_DBL_MIN  % .16e\n", sqrt (DBL_MIN));
  printf ("#define GSL_ROOT3_DBL_MIN % .16e\n", rt3 (DBL_MIN));
  printf ("#define GSL_ROOT4_DBL_MIN % .16e\n", rt4 (DBL_MIN));
  printf ("#define GSL_ROOT5_DBL_MIN % .16e\n", rt5 (DBL_MIN));
  printf ("#define GSL_ROOT6_DBL_MIN % .16e\n", rt6 (DBL_MIN));
  printf ("#define GSL_LOG_DBL_MIN   % .16e\n", log (DBL_MIN));
  printf ("\n");

  printf ("#define GSL_DBL_MAX       % .16e\n", DBL_MAX);
  printf ("#define GSL_SQRT_DBL_MAX  % .16e\n", sqrt (DBL_MAX));
  printf ("#define GSL_ROOT3_DBL_MAX % .16e\n", rt3 (DBL_MAX));
  printf ("#define GSL_ROOT4_DBL_MAX % .16e\n", rt4 (DBL_MAX));
  printf ("#define GSL_ROOT5_DBL_MAX % .16e\n", rt5 (DBL_MAX));
  printf ("#define GSL_ROOT6_DBL_MAX % .16e\n", rt6 (DBL_MAX));
  printf ("#define GSL_LOG_DBL_MAX   % .16e\n", log (DBL_MAX));
  printf ("\n");

  printf ("#define GSL_FLT_EPSILON       % .16e\n", FLT_EPSILON);
  printf ("#define GSL_SQRT_FLT_EPSILON  % .16e\n", sqrt (FLT_EPSILON));
  printf ("#define GSL_ROOT3_FLT_EPSILON % .16e\n", rt3 (FLT_EPSILON));
  printf ("#define GSL_ROOT4_FLT_EPSILON % .16e\n", rt4 (FLT_EPSILON));
  printf ("#define GSL_ROOT5_FLT_EPSILON % .16e\n", rt5 (FLT_EPSILON));
  printf ("#define GSL_ROOT6_FLT_EPSILON % .16e\n", rt6 (FLT_EPSILON));
  printf ("#define GSL_LOG_FLT_EPSILON   % .16e\n", log (FLT_EPSILON));
  printf ("\n");

  printf ("#define GSL_FLT_MIN       % .16e\n", FLT_MIN);
  printf ("#define GSL_SQRT_FLT_MIN  % .16e\n", sqrt (FLT_MIN));
  printf ("#define GSL_ROOT3_FLT_MIN % .16e\n", rt3 (FLT_MIN));
  printf ("#define GSL_ROOT4_FLT_MIN % .16e\n", rt4 (FLT_MIN));
  printf ("#define GSL_ROOT5_FLT_MIN % .16e\n", rt5 (FLT_MIN));
  printf ("#define GSL_ROOT6_FLT_MIN % .16e\n", rt6 (FLT_MIN));
  printf ("#define GSL_LOG_FLT_MIN   % .16e\n", log (FLT_MIN));
  printf ("\n");

  printf ("#define GSL_FLT_MAX       % .16e\n", FLT_MAX);
  printf ("#define GSL_SQRT_FLT_MAX  % .16e\n", sqrt (FLT_MAX));
  printf ("#define GSL_ROOT3_FLT_MAX % .16e\n", rt3 (FLT_MAX));
  printf ("#define GSL_ROOT4_FLT_MAX % .16e\n", rt4 (FLT_MAX));
  printf ("#define GSL_ROOT5_FLT_MAX % .16e\n", rt5 (FLT_MAX));
  printf ("#define GSL_ROOT6_FLT_MAX % .16e\n", rt6 (FLT_MAX));
  printf ("#define GSL_LOG_FLT_MAX   % .16e\n", log (FLT_MAX));
  printf ("\n");

  return 0;
}
