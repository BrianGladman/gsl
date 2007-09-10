/* sys/params.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */
#include <config.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#define rt3(x) pow((x), 1.0 / 3.0)
#define rt4(x) pow((x), 1.0 / 4.0)
#define rt5(x) pow((x), 1.0 / 5.0)
#define rt6(x) pow((x), 1.0 / 6.0)

#define SFLT_EPSILON 0.00048828125  /* 2^(-11) */


int
main (void)
{
  printf ("#define GSL_DBL_EPSILON       % .16e\n", DBL_EPSILON);
  printf ("#define GSL_SQRT_DBL_EPSILON  % .16e\n", sqrt (DBL_EPSILON));
  printf ("#define GSL_ROOT3_DBL_EPSILON % .16e\n", rt3 (DBL_EPSILON));
  printf ("#define GSL_ROOT4_DBL_EPSILON % .16e\n", rt4 (DBL_EPSILON));
  printf ("#define GSL_ROOT5_DBL_EPSILON % .16e\n", rt5 (DBL_EPSILON));
  printf ("#define GSL_ROOT6_DBL_EPSILON % .16e\n", rt6 (DBL_EPSILON));
  printf ("#define GSL_LOG_DBL_EPSILON  (% .16e)\n", log (DBL_EPSILON));
  printf ("\n");

  printf ("#define GSL_DBL_MIN       % .16e\n", DBL_MIN);
  printf ("#define GSL_SQRT_DBL_MIN  % .16e\n", sqrt (DBL_MIN));
  printf ("#define GSL_ROOT3_DBL_MIN % .16e\n", rt3 (DBL_MIN));
  printf ("#define GSL_ROOT4_DBL_MIN % .16e\n", rt4 (DBL_MIN));
  printf ("#define GSL_ROOT5_DBL_MIN % .16e\n", rt5 (DBL_MIN));
  printf ("#define GSL_ROOT6_DBL_MIN % .16e\n", rt6 (DBL_MIN));
  printf ("#define GSL_LOG_DBL_MIN  (% .16e)\n", log (DBL_MIN));
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
  printf ("#define GSL_LOG_FLT_EPSILON  (% .16e)\n", log (FLT_EPSILON));
  printf ("\n");

  printf ("#define GSL_FLT_MIN       % .16e\n", FLT_MIN);
  printf ("#define GSL_SQRT_FLT_MIN  % .16e\n", sqrt (FLT_MIN));
  printf ("#define GSL_ROOT3_FLT_MIN % .16e\n", rt3 (FLT_MIN));
  printf ("#define GSL_ROOT4_FLT_MIN % .16e\n", rt4 (FLT_MIN));
  printf ("#define GSL_ROOT5_FLT_MIN % .16e\n", rt5 (FLT_MIN));
  printf ("#define GSL_ROOT6_FLT_MIN % .16e\n", rt6 (FLT_MIN));
  printf ("#define GSL_LOG_FLT_MIN  (% .16e)\n", log (FLT_MIN));
  printf ("\n");

  printf ("#define GSL_FLT_MAX       % .16e\n", FLT_MAX);
  printf ("#define GSL_SQRT_FLT_MAX  % .16e\n", sqrt (FLT_MAX));
  printf ("#define GSL_ROOT3_FLT_MAX % .16e\n", rt3 (FLT_MAX));
  printf ("#define GSL_ROOT4_FLT_MAX % .16e\n", rt4 (FLT_MAX));
  printf ("#define GSL_ROOT5_FLT_MAX % .16e\n", rt5 (FLT_MAX));
  printf ("#define GSL_ROOT6_FLT_MAX % .16e\n", rt6 (FLT_MAX));
  printf ("#define GSL_LOG_FLT_MAX   % .16e\n", log (FLT_MAX));
  printf ("\n");

  printf ("#define GSL_SFLT_EPSILON       % .16e\n", SFLT_EPSILON);
  printf ("#define GSL_SQRT_SFLT_EPSILON  % .16e\n", sqrt (SFLT_EPSILON));
  printf ("#define GSL_ROOT3_SFLT_EPSILON % .16e\n", rt3 (SFLT_EPSILON));
  printf ("#define GSL_ROOT4_SFLT_EPSILON % .16e\n", rt4 (SFLT_EPSILON));
  printf ("#define GSL_ROOT5_SFLT_EPSILON % .16e\n", rt5 (SFLT_EPSILON));
  printf ("#define GSL_ROOT6_SFLT_EPSILON % .16e\n", rt6 (SFLT_EPSILON));
  printf ("#define GSL_LOG_SFLT_EPSILON  (% .16e)\n", log (SFLT_EPSILON));
  printf ("\n");

  return 0;
}

/* This is the output on a PA-RISC HPUX-9 system

#define GSL_DBL_EPSILON        2.2204460492503131e-16
#define GSL_SQRT_DBL_EPSILON   1.4901161193847656e-08
#define GSL_ROOT3_DBL_EPSILON  6.0554544523933429e-06
#define GSL_ROOT4_DBL_EPSILON  1.2207031250000000e-04
#define GSL_ROOT5_DBL_EPSILON  7.4009597974140505e-04
#define GSL_ROOT6_DBL_EPSILON  2.4607833005759251e-03
#define GSL_LOG_DBL_EPSILON   -3.6043653389117154e+01

#define GSL_DBL_MIN        2.2250738585072014e-308
#define GSL_SQRT_DBL_MIN   1.4916681462400413e-154
#define GSL_ROOT3_DBL_MIN  2.8126442852362996e-103
#define GSL_ROOT4_DBL_MIN  1.2213386697554620e-77
#define GSL_ROOT5_DBL_MIN  2.9476022969691763e-62
#define GSL_ROOT6_DBL_MIN  5.3034368905798218e-52
#define GSL_LOG_DBL_MIN   -7.0839641853226408e+02

#define GSL_DBL_MAX        1.7976931348623157e+308
#define GSL_SQRT_DBL_MAX   1.3407807929942596e+154
#define GSL_ROOT3_DBL_MAX  5.6438030941222897e+102
#define GSL_ROOT4_DBL_MAX  1.1579208923731620e+77
#define GSL_ROOT5_DBL_MAX  4.4765466227572707e+61
#define GSL_ROOT6_DBL_MAX  2.3756689782295612e+51
#define GSL_LOG_DBL_MAX    7.0978271289338397e+02

#define GSL_FLT_EPSILON        1.1920928955078125e-07
#define GSL_SQRT_FLT_EPSILON   3.4526698300124393e-04
#define GSL_ROOT3_FLT_EPSILON  4.9215666011518501e-03
#define GSL_ROOT4_FLT_EPSILON  1.8581361171917516e-02
#define GSL_ROOT5_FLT_EPSILON  4.1234622211652937e-02
#define GSL_ROOT6_FLT_EPSILON  7.0153878019335827e-02
#define GSL_LOG_FLT_EPSILON   -1.5942385152878742e+01

#define GSL_FLT_MIN        1.1754943508222875e-38
#define GSL_SQRT_FLT_MIN   1.0842021724855044e-19
#define GSL_ROOT3_FLT_MIN  2.2737367544323241e-13
#define GSL_ROOT4_FLT_MIN  3.2927225399135965e-10
#define GSL_ROOT5_FLT_MIN  2.5944428542140822e-08
#define GSL_ROOT6_FLT_MIN  4.7683715820312542e-07
#define GSL_LOG_FLT_MIN   -8.7336544750553102e+01

#define GSL_FLT_MAX        3.4028234663852886e+38
#define GSL_SQRT_FLT_MAX   1.8446743523953730e+19
#define GSL_ROOT3_FLT_MAX  6.9814635196223242e+12
#define GSL_ROOT4_FLT_MAX  4.2949672319999986e+09
#define GSL_ROOT5_FLT_MAX  5.0859007855960041e+07
#define GSL_ROOT6_FLT_MAX  2.6422459233807749e+06
#define GSL_LOG_FLT_MAX    8.8722839052068352e+01

*/



