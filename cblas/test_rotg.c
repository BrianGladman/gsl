#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

void
test_rotg () {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   float a = -1.5;
   float b = -1.5;
   float c;
   float s;
   float r_expected = 2.12132034356;
   float z_expected = -1.41421356237;
   float c_expected = -0.707106781187;
   float s_expected = -0.707106781187;
   cblas_srotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, flteps, "srotg(case 56)");
   gsl_test_rel(b, z_expected, flteps, "srotg(case 57)");
   gsl_test_rel(c, c_expected, flteps, "srotg(case 58)");
   gsl_test_rel(s, s_expected, flteps, "srotg(case 59)");
  };


  {
   float a = -1.5;
   float b = -1;
   float c;
   float s;
   float r_expected = -1.80277563773;
   float z_expected = 0.554700196225;
   float c_expected = 0.832050294338;
   float s_expected = 0.554700196225;
   cblas_srotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, flteps, "srotg(case 60)");
   gsl_test_rel(b, z_expected, flteps, "srotg(case 61)");
   gsl_test_rel(c, c_expected, flteps, "srotg(case 62)");
   gsl_test_rel(s, s_expected, flteps, "srotg(case 63)");
  };


  {
   float a = -1.5;
   float b = -0.1;
   float c;
   float s;
   float r_expected = -1.50332963784;
   float z_expected = 0.0665190105238;
   float c_expected = 0.997785157857;
   float s_expected = 0.0665190105238;
   cblas_srotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, flteps, "srotg(case 64)");
   gsl_test_rel(b, z_expected, flteps, "srotg(case 65)");
   gsl_test_rel(c, c_expected, flteps, "srotg(case 66)");
   gsl_test_rel(s, s_expected, flteps, "srotg(case 67)");
  };


  {
   float a = -1.5;
   float b = 0;
   float c;
   float s;
   float r_expected = -1.5;
   float z_expected = 0;
   float c_expected = 1;
   float s_expected = 0;
   cblas_srotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, flteps, "srotg(case 68)");
   gsl_test_rel(b, z_expected, flteps, "srotg(case 69)");
   gsl_test_rel(c, c_expected, flteps, "srotg(case 70)");
   gsl_test_rel(s, s_expected, flteps, "srotg(case 71)");
  };


  {
   float a = -1.5;
   float b = 0.1;
   float c;
   float s;
   float r_expected = -1.50332963784;
   float z_expected = -0.0665190105238;
   float c_expected = 0.997785157857;
   float s_expected = -0.0665190105238;
   cblas_srotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, flteps, "srotg(case 72)");
   gsl_test_rel(b, z_expected, flteps, "srotg(case 73)");
   gsl_test_rel(c, c_expected, flteps, "srotg(case 74)");
   gsl_test_rel(s, s_expected, flteps, "srotg(case 75)");
  };


  {
   float a = -1.5;
   float b = 1;
   float c;
   float s;
   float r_expected = -1.80277563773;
   float z_expected = -0.554700196225;
   float c_expected = 0.832050294338;
   float s_expected = -0.554700196225;
   cblas_srotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, flteps, "srotg(case 76)");
   gsl_test_rel(b, z_expected, flteps, "srotg(case 77)");
   gsl_test_rel(c, c_expected, flteps, "srotg(case 78)");
   gsl_test_rel(s, s_expected, flteps, "srotg(case 79)");
  };


  {
   float a = -1.5;
   float b = 1.5;
   float c;
   float s;
   float r_expected = 2.12132034356;
   float z_expected = -1.41421356237;
   float c_expected = -0.707106781187;
   float s_expected = 0.707106781187;
   cblas_srotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, flteps, "srotg(case 80)");
   gsl_test_rel(b, z_expected, flteps, "srotg(case 81)");
   gsl_test_rel(c, c_expected, flteps, "srotg(case 82)");
   gsl_test_rel(s, s_expected, flteps, "srotg(case 83)");
  };


  {
   float a = -1;
   float b = -1.5;
   float c;
   float s;
   float r_expected = 1.80277563773;
   float z_expected = -1.80277563773;
   float c_expected = -0.554700196225;
   float s_expected = -0.832050294338;
   cblas_srotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, flteps, "srotg(case 84)");
   gsl_test_rel(b, z_expected, flteps, "srotg(case 85)");
   gsl_test_rel(c, c_expected, flteps, "srotg(case 86)");
   gsl_test_rel(s, s_expected, flteps, "srotg(case 87)");
  };


  {
   float a = -1;
   float b = -1;
   float c;
   float s;
   float r_expected = 1.41421356237;
   float z_expected = -1.41421356237;
   float c_expected = -0.707106781187;
   float s_expected = -0.707106781187;
   cblas_srotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, flteps, "srotg(case 88)");
   gsl_test_rel(b, z_expected, flteps, "srotg(case 89)");
   gsl_test_rel(c, c_expected, flteps, "srotg(case 90)");
   gsl_test_rel(s, s_expected, flteps, "srotg(case 91)");
  };


  {
   float a = -1;
   float b = -0.1;
   float c;
   float s;
   float r_expected = -1.00498756211;
   float z_expected = 0.099503719021;
   float c_expected = 0.99503719021;
   float s_expected = 0.099503719021;
   cblas_srotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, flteps, "srotg(case 92)");
   gsl_test_rel(b, z_expected, flteps, "srotg(case 93)");
   gsl_test_rel(c, c_expected, flteps, "srotg(case 94)");
   gsl_test_rel(s, s_expected, flteps, "srotg(case 95)");
  };


  {
   float a = -1;
   float b = 0;
   float c;
   float s;
   float r_expected = -1;
   float z_expected = 0;
   float c_expected = 1;
   float s_expected = 0;
   cblas_srotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, flteps, "srotg(case 96)");
   gsl_test_rel(b, z_expected, flteps, "srotg(case 97)");
   gsl_test_rel(c, c_expected, flteps, "srotg(case 98)");
   gsl_test_rel(s, s_expected, flteps, "srotg(case 99)");
  };


  {
   float a = -1;
   float b = 0.1;
   float c;
   float s;
   float r_expected = -1.00498756211;
   float z_expected = -0.099503719021;
   float c_expected = 0.99503719021;
   float s_expected = -0.099503719021;
   cblas_srotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, flteps, "srotg(case 100)");
   gsl_test_rel(b, z_expected, flteps, "srotg(case 101)");
   gsl_test_rel(c, c_expected, flteps, "srotg(case 102)");
   gsl_test_rel(s, s_expected, flteps, "srotg(case 103)");
  };


  {
   float a = -1;
   float b = 1;
   float c;
   float s;
   float r_expected = 1.41421356237;
   float z_expected = -1.41421356237;
   float c_expected = -0.707106781187;
   float s_expected = 0.707106781187;
   cblas_srotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, flteps, "srotg(case 104)");
   gsl_test_rel(b, z_expected, flteps, "srotg(case 105)");
   gsl_test_rel(c, c_expected, flteps, "srotg(case 106)");
   gsl_test_rel(s, s_expected, flteps, "srotg(case 107)");
  };


  {
   float a = -1;
   float b = 1.5;
   float c;
   float s;
   float r_expected = 1.80277563773;
   float z_expected = -1.80277563773;
   float c_expected = -0.554700196225;
   float s_expected = 0.832050294338;
   cblas_srotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, flteps, "srotg(case 108)");
   gsl_test_rel(b, z_expected, flteps, "srotg(case 109)");
   gsl_test_rel(c, c_expected, flteps, "srotg(case 110)");
   gsl_test_rel(s, s_expected, flteps, "srotg(case 111)");
  };


  {
   float a = -0.1;
   float b = -1.5;
   float c;
   float s;
   float r_expected = 1.50332963784;
   float z_expected = -15.0332963784;
   float c_expected = -0.0665190105238;
   float s_expected = -0.997785157857;
   cblas_srotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, flteps, "srotg(case 112)");
   gsl_test_rel(b, z_expected, flteps, "srotg(case 113)");
   gsl_test_rel(c, c_expected, flteps, "srotg(case 114)");
   gsl_test_rel(s, s_expected, flteps, "srotg(case 115)");
  };


  {
   float a = -0.1;
   float b = -1;
   float c;
   float s;
   float r_expected = 1.00498756211;
   float z_expected = -10.0498756211;
   float c_expected = -0.099503719021;
   float s_expected = -0.99503719021;
   cblas_srotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, flteps, "srotg(case 116)");
   gsl_test_rel(b, z_expected, flteps, "srotg(case 117)");
   gsl_test_rel(c, c_expected, flteps, "srotg(case 118)");
   gsl_test_rel(s, s_expected, flteps, "srotg(case 119)");
  };


  {
   float a = -0.1;
   float b = -0.1;
   float c;
   float s;
   float r_expected = 0.141421356237;
   float z_expected = -1.41421356237;
   float c_expected = -0.707106781187;
   float s_expected = -0.707106781187;
   cblas_srotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, flteps, "srotg(case 120)");
   gsl_test_rel(b, z_expected, flteps, "srotg(case 121)");
   gsl_test_rel(c, c_expected, flteps, "srotg(case 122)");
   gsl_test_rel(s, s_expected, flteps, "srotg(case 123)");
  };


  {
   float a = -0.1;
   float b = 0;
   float c;
   float s;
   float r_expected = -0.1;
   float z_expected = 0;
   float c_expected = 1;
   float s_expected = 0;
   cblas_srotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, flteps, "srotg(case 124)");
   gsl_test_rel(b, z_expected, flteps, "srotg(case 125)");
   gsl_test_rel(c, c_expected, flteps, "srotg(case 126)");
   gsl_test_rel(s, s_expected, flteps, "srotg(case 127)");
  };


  {
   float a = -0.1;
   float b = 0.1;
   float c;
   float s;
   float r_expected = 0.141421356237;
   float z_expected = -1.41421356237;
   float c_expected = -0.707106781187;
   float s_expected = 0.707106781187;
   cblas_srotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, flteps, "srotg(case 128)");
   gsl_test_rel(b, z_expected, flteps, "srotg(case 129)");
   gsl_test_rel(c, c_expected, flteps, "srotg(case 130)");
   gsl_test_rel(s, s_expected, flteps, "srotg(case 131)");
  };


  {
   float a = -0.1;
   float b = 1;
   float c;
   float s;
   float r_expected = 1.00498756211;
   float z_expected = -10.0498756211;
   float c_expected = -0.099503719021;
   float s_expected = 0.99503719021;
   cblas_srotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, flteps, "srotg(case 132)");
   gsl_test_rel(b, z_expected, flteps, "srotg(case 133)");
   gsl_test_rel(c, c_expected, flteps, "srotg(case 134)");
   gsl_test_rel(s, s_expected, flteps, "srotg(case 135)");
  };


  {
   float a = -0.1;
   float b = 1.5;
   float c;
   float s;
   float r_expected = 1.50332963784;
   float z_expected = -15.0332963784;
   float c_expected = -0.0665190105238;
   float s_expected = 0.997785157857;
   cblas_srotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, flteps, "srotg(case 136)");
   gsl_test_rel(b, z_expected, flteps, "srotg(case 137)");
   gsl_test_rel(c, c_expected, flteps, "srotg(case 138)");
   gsl_test_rel(s, s_expected, flteps, "srotg(case 139)");
  };


  {
   float a = 0;
   float b = -1.5;
   float c;
   float s;
   float r_expected = -1.5;
   float z_expected = 1;
   float c_expected = 0;
   float s_expected = 1;
   cblas_srotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, flteps, "srotg(case 140)");
   gsl_test_rel(b, z_expected, flteps, "srotg(case 141)");
   gsl_test_rel(c, c_expected, flteps, "srotg(case 142)");
   gsl_test_rel(s, s_expected, flteps, "srotg(case 143)");
  };


  {
   float a = 0;
   float b = -1;
   float c;
   float s;
   float r_expected = -1;
   float z_expected = 1;
   float c_expected = 0;
   float s_expected = 1;
   cblas_srotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, flteps, "srotg(case 144)");
   gsl_test_rel(b, z_expected, flteps, "srotg(case 145)");
   gsl_test_rel(c, c_expected, flteps, "srotg(case 146)");
   gsl_test_rel(s, s_expected, flteps, "srotg(case 147)");
  };


  {
   float a = 0;
   float b = -0.1;
   float c;
   float s;
   float r_expected = -0.1;
   float z_expected = 1;
   float c_expected = 0;
   float s_expected = 1;
   cblas_srotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, flteps, "srotg(case 148)");
   gsl_test_rel(b, z_expected, flteps, "srotg(case 149)");
   gsl_test_rel(c, c_expected, flteps, "srotg(case 150)");
   gsl_test_rel(s, s_expected, flteps, "srotg(case 151)");
  };


  {
   float a = 0;
   float b = 0;
   float c;
   float s;
   float r_expected = 0;
   float z_expected = 1;
   float c_expected = 1;
   float s_expected = 0;
   cblas_srotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, flteps, "srotg(case 152)");
   gsl_test_rel(b, z_expected, flteps, "srotg(case 153)");
   gsl_test_rel(c, c_expected, flteps, "srotg(case 154)");
   gsl_test_rel(s, s_expected, flteps, "srotg(case 155)");
  };


  {
   float a = 0;
   float b = 0.1;
   float c;
   float s;
   float r_expected = 0.1;
   float z_expected = 1;
   float c_expected = 0;
   float s_expected = 1;
   cblas_srotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, flteps, "srotg(case 156)");
   gsl_test_rel(b, z_expected, flteps, "srotg(case 157)");
   gsl_test_rel(c, c_expected, flteps, "srotg(case 158)");
   gsl_test_rel(s, s_expected, flteps, "srotg(case 159)");
  };


  {
   float a = 0;
   float b = 1;
   float c;
   float s;
   float r_expected = 1;
   float z_expected = 1;
   float c_expected = 0;
   float s_expected = 1;
   cblas_srotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, flteps, "srotg(case 160)");
   gsl_test_rel(b, z_expected, flteps, "srotg(case 161)");
   gsl_test_rel(c, c_expected, flteps, "srotg(case 162)");
   gsl_test_rel(s, s_expected, flteps, "srotg(case 163)");
  };


  {
   float a = 0;
   float b = 1.5;
   float c;
   float s;
   float r_expected = 1.5;
   float z_expected = 1;
   float c_expected = 0;
   float s_expected = 1;
   cblas_srotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, flteps, "srotg(case 164)");
   gsl_test_rel(b, z_expected, flteps, "srotg(case 165)");
   gsl_test_rel(c, c_expected, flteps, "srotg(case 166)");
   gsl_test_rel(s, s_expected, flteps, "srotg(case 167)");
  };


  {
   float a = 0.1;
   float b = -1.5;
   float c;
   float s;
   float r_expected = 1.50332963784;
   float z_expected = 15.0332963784;
   float c_expected = 0.0665190105238;
   float s_expected = -0.997785157857;
   cblas_srotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, flteps, "srotg(case 168)");
   gsl_test_rel(b, z_expected, flteps, "srotg(case 169)");
   gsl_test_rel(c, c_expected, flteps, "srotg(case 170)");
   gsl_test_rel(s, s_expected, flteps, "srotg(case 171)");
  };


  {
   float a = 0.1;
   float b = -1;
   float c;
   float s;
   float r_expected = 1.00498756211;
   float z_expected = 10.0498756211;
   float c_expected = 0.099503719021;
   float s_expected = -0.99503719021;
   cblas_srotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, flteps, "srotg(case 172)");
   gsl_test_rel(b, z_expected, flteps, "srotg(case 173)");
   gsl_test_rel(c, c_expected, flteps, "srotg(case 174)");
   gsl_test_rel(s, s_expected, flteps, "srotg(case 175)");
  };


  {
   float a = 0.1;
   float b = -0.1;
   float c;
   float s;
   float r_expected = 0.141421356237;
   float z_expected = 1.41421356237;
   float c_expected = 0.707106781187;
   float s_expected = -0.707106781187;
   cblas_srotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, flteps, "srotg(case 176)");
   gsl_test_rel(b, z_expected, flteps, "srotg(case 177)");
   gsl_test_rel(c, c_expected, flteps, "srotg(case 178)");
   gsl_test_rel(s, s_expected, flteps, "srotg(case 179)");
  };


  {
   float a = 0.1;
   float b = 0;
   float c;
   float s;
   float r_expected = 0.1;
   float z_expected = 0;
   float c_expected = 1;
   float s_expected = 0;
   cblas_srotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, flteps, "srotg(case 180)");
   gsl_test_rel(b, z_expected, flteps, "srotg(case 181)");
   gsl_test_rel(c, c_expected, flteps, "srotg(case 182)");
   gsl_test_rel(s, s_expected, flteps, "srotg(case 183)");
  };


  {
   float a = 0.1;
   float b = 0.1;
   float c;
   float s;
   float r_expected = 0.141421356237;
   float z_expected = 1.41421356237;
   float c_expected = 0.707106781187;
   float s_expected = 0.707106781187;
   cblas_srotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, flteps, "srotg(case 184)");
   gsl_test_rel(b, z_expected, flteps, "srotg(case 185)");
   gsl_test_rel(c, c_expected, flteps, "srotg(case 186)");
   gsl_test_rel(s, s_expected, flteps, "srotg(case 187)");
  };


  {
   float a = 0.1;
   float b = 1;
   float c;
   float s;
   float r_expected = 1.00498756211;
   float z_expected = 10.0498756211;
   float c_expected = 0.099503719021;
   float s_expected = 0.99503719021;
   cblas_srotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, flteps, "srotg(case 188)");
   gsl_test_rel(b, z_expected, flteps, "srotg(case 189)");
   gsl_test_rel(c, c_expected, flteps, "srotg(case 190)");
   gsl_test_rel(s, s_expected, flteps, "srotg(case 191)");
  };


  {
   float a = 0.1;
   float b = 1.5;
   float c;
   float s;
   float r_expected = 1.50332963784;
   float z_expected = 15.0332963784;
   float c_expected = 0.0665190105238;
   float s_expected = 0.997785157857;
   cblas_srotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, flteps, "srotg(case 192)");
   gsl_test_rel(b, z_expected, flteps, "srotg(case 193)");
   gsl_test_rel(c, c_expected, flteps, "srotg(case 194)");
   gsl_test_rel(s, s_expected, flteps, "srotg(case 195)");
  };


  {
   float a = 1;
   float b = -1.5;
   float c;
   float s;
   float r_expected = 1.80277563773;
   float z_expected = 1.80277563773;
   float c_expected = 0.554700196225;
   float s_expected = -0.832050294338;
   cblas_srotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, flteps, "srotg(case 196)");
   gsl_test_rel(b, z_expected, flteps, "srotg(case 197)");
   gsl_test_rel(c, c_expected, flteps, "srotg(case 198)");
   gsl_test_rel(s, s_expected, flteps, "srotg(case 199)");
  };


  {
   float a = 1;
   float b = -1;
   float c;
   float s;
   float r_expected = 1.41421356237;
   float z_expected = 1.41421356237;
   float c_expected = 0.707106781187;
   float s_expected = -0.707106781187;
   cblas_srotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, flteps, "srotg(case 200)");
   gsl_test_rel(b, z_expected, flteps, "srotg(case 201)");
   gsl_test_rel(c, c_expected, flteps, "srotg(case 202)");
   gsl_test_rel(s, s_expected, flteps, "srotg(case 203)");
  };


  {
   float a = 1;
   float b = -0.1;
   float c;
   float s;
   float r_expected = 1.00498756211;
   float z_expected = -0.099503719021;
   float c_expected = 0.99503719021;
   float s_expected = -0.099503719021;
   cblas_srotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, flteps, "srotg(case 204)");
   gsl_test_rel(b, z_expected, flteps, "srotg(case 205)");
   gsl_test_rel(c, c_expected, flteps, "srotg(case 206)");
   gsl_test_rel(s, s_expected, flteps, "srotg(case 207)");
  };


  {
   float a = 1;
   float b = 0;
   float c;
   float s;
   float r_expected = 1;
   float z_expected = 0;
   float c_expected = 1;
   float s_expected = 0;
   cblas_srotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, flteps, "srotg(case 208)");
   gsl_test_rel(b, z_expected, flteps, "srotg(case 209)");
   gsl_test_rel(c, c_expected, flteps, "srotg(case 210)");
   gsl_test_rel(s, s_expected, flteps, "srotg(case 211)");
  };


  {
   float a = 1;
   float b = 0.1;
   float c;
   float s;
   float r_expected = 1.00498756211;
   float z_expected = 0.099503719021;
   float c_expected = 0.99503719021;
   float s_expected = 0.099503719021;
   cblas_srotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, flteps, "srotg(case 212)");
   gsl_test_rel(b, z_expected, flteps, "srotg(case 213)");
   gsl_test_rel(c, c_expected, flteps, "srotg(case 214)");
   gsl_test_rel(s, s_expected, flteps, "srotg(case 215)");
  };


  {
   float a = 1;
   float b = 1;
   float c;
   float s;
   float r_expected = 1.41421356237;
   float z_expected = 1.41421356237;
   float c_expected = 0.707106781187;
   float s_expected = 0.707106781187;
   cblas_srotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, flteps, "srotg(case 216)");
   gsl_test_rel(b, z_expected, flteps, "srotg(case 217)");
   gsl_test_rel(c, c_expected, flteps, "srotg(case 218)");
   gsl_test_rel(s, s_expected, flteps, "srotg(case 219)");
  };


  {
   float a = 1;
   float b = 1.5;
   float c;
   float s;
   float r_expected = 1.80277563773;
   float z_expected = 1.80277563773;
   float c_expected = 0.554700196225;
   float s_expected = 0.832050294338;
   cblas_srotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, flteps, "srotg(case 220)");
   gsl_test_rel(b, z_expected, flteps, "srotg(case 221)");
   gsl_test_rel(c, c_expected, flteps, "srotg(case 222)");
   gsl_test_rel(s, s_expected, flteps, "srotg(case 223)");
  };


  {
   float a = 1.5;
   float b = -1.5;
   float c;
   float s;
   float r_expected = 2.12132034356;
   float z_expected = 1.41421356237;
   float c_expected = 0.707106781187;
   float s_expected = -0.707106781187;
   cblas_srotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, flteps, "srotg(case 224)");
   gsl_test_rel(b, z_expected, flteps, "srotg(case 225)");
   gsl_test_rel(c, c_expected, flteps, "srotg(case 226)");
   gsl_test_rel(s, s_expected, flteps, "srotg(case 227)");
  };


  {
   float a = 1.5;
   float b = -1;
   float c;
   float s;
   float r_expected = 1.80277563773;
   float z_expected = -0.554700196225;
   float c_expected = 0.832050294338;
   float s_expected = -0.554700196225;
   cblas_srotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, flteps, "srotg(case 228)");
   gsl_test_rel(b, z_expected, flteps, "srotg(case 229)");
   gsl_test_rel(c, c_expected, flteps, "srotg(case 230)");
   gsl_test_rel(s, s_expected, flteps, "srotg(case 231)");
  };


  {
   float a = 1.5;
   float b = -0.1;
   float c;
   float s;
   float r_expected = 1.50332963784;
   float z_expected = -0.0665190105238;
   float c_expected = 0.997785157857;
   float s_expected = -0.0665190105238;
   cblas_srotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, flteps, "srotg(case 232)");
   gsl_test_rel(b, z_expected, flteps, "srotg(case 233)");
   gsl_test_rel(c, c_expected, flteps, "srotg(case 234)");
   gsl_test_rel(s, s_expected, flteps, "srotg(case 235)");
  };


  {
   float a = 1.5;
   float b = 0;
   float c;
   float s;
   float r_expected = 1.5;
   float z_expected = 0;
   float c_expected = 1;
   float s_expected = 0;
   cblas_srotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, flteps, "srotg(case 236)");
   gsl_test_rel(b, z_expected, flteps, "srotg(case 237)");
   gsl_test_rel(c, c_expected, flteps, "srotg(case 238)");
   gsl_test_rel(s, s_expected, flteps, "srotg(case 239)");
  };


  {
   float a = 1.5;
   float b = 0.1;
   float c;
   float s;
   float r_expected = 1.50332963784;
   float z_expected = 0.0665190105238;
   float c_expected = 0.997785157857;
   float s_expected = 0.0665190105238;
   cblas_srotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, flteps, "srotg(case 240)");
   gsl_test_rel(b, z_expected, flteps, "srotg(case 241)");
   gsl_test_rel(c, c_expected, flteps, "srotg(case 242)");
   gsl_test_rel(s, s_expected, flteps, "srotg(case 243)");
  };


  {
   float a = 1.5;
   float b = 1;
   float c;
   float s;
   float r_expected = 1.80277563773;
   float z_expected = 0.554700196225;
   float c_expected = 0.832050294338;
   float s_expected = 0.554700196225;
   cblas_srotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, flteps, "srotg(case 244)");
   gsl_test_rel(b, z_expected, flteps, "srotg(case 245)");
   gsl_test_rel(c, c_expected, flteps, "srotg(case 246)");
   gsl_test_rel(s, s_expected, flteps, "srotg(case 247)");
  };


  {
   float a = 1.5;
   float b = 1.5;
   float c;
   float s;
   float r_expected = 2.12132034356;
   float z_expected = 1.41421356237;
   float c_expected = 0.707106781187;
   float s_expected = 0.707106781187;
   cblas_srotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, flteps, "srotg(case 248)");
   gsl_test_rel(b, z_expected, flteps, "srotg(case 249)");
   gsl_test_rel(c, c_expected, flteps, "srotg(case 250)");
   gsl_test_rel(s, s_expected, flteps, "srotg(case 251)");
  };


  {
   double a = -1.5;
   double b = -1.5;
   double c;
   double s;
   double r_expected = 2.12132034356;
   double z_expected = -1.41421356237;
   double c_expected = -0.707106781187;
   double s_expected = -0.707106781187;
   cblas_drotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, dbleps, "drotg(case 252)");
   gsl_test_rel(b, z_expected, dbleps, "drotg(case 253)");
   gsl_test_rel(c, c_expected, dbleps, "drotg(case 254)");
   gsl_test_rel(s, s_expected, dbleps, "drotg(case 255)");
  };


  {
   double a = -1.5;
   double b = -1;
   double c;
   double s;
   double r_expected = -1.80277563773;
   double z_expected = 0.554700196225;
   double c_expected = 0.832050294338;
   double s_expected = 0.554700196225;
   cblas_drotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, dbleps, "drotg(case 256)");
   gsl_test_rel(b, z_expected, dbleps, "drotg(case 257)");
   gsl_test_rel(c, c_expected, dbleps, "drotg(case 258)");
   gsl_test_rel(s, s_expected, dbleps, "drotg(case 259)");
  };


  {
   double a = -1.5;
   double b = -0.1;
   double c;
   double s;
   double r_expected = -1.50332963784;
   double z_expected = 0.0665190105238;
   double c_expected = 0.997785157857;
   double s_expected = 0.0665190105238;
   cblas_drotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, dbleps, "drotg(case 260)");
   gsl_test_rel(b, z_expected, dbleps, "drotg(case 261)");
   gsl_test_rel(c, c_expected, dbleps, "drotg(case 262)");
   gsl_test_rel(s, s_expected, dbleps, "drotg(case 263)");
  };


  {
   double a = -1.5;
   double b = 0;
   double c;
   double s;
   double r_expected = -1.5;
   double z_expected = 0;
   double c_expected = 1;
   double s_expected = 0;
   cblas_drotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, dbleps, "drotg(case 264)");
   gsl_test_rel(b, z_expected, dbleps, "drotg(case 265)");
   gsl_test_rel(c, c_expected, dbleps, "drotg(case 266)");
   gsl_test_rel(s, s_expected, dbleps, "drotg(case 267)");
  };


  {
   double a = -1.5;
   double b = 0.1;
   double c;
   double s;
   double r_expected = -1.50332963784;
   double z_expected = -0.0665190105238;
   double c_expected = 0.997785157857;
   double s_expected = -0.0665190105238;
   cblas_drotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, dbleps, "drotg(case 268)");
   gsl_test_rel(b, z_expected, dbleps, "drotg(case 269)");
   gsl_test_rel(c, c_expected, dbleps, "drotg(case 270)");
   gsl_test_rel(s, s_expected, dbleps, "drotg(case 271)");
  };


  {
   double a = -1.5;
   double b = 1;
   double c;
   double s;
   double r_expected = -1.80277563773;
   double z_expected = -0.554700196225;
   double c_expected = 0.832050294338;
   double s_expected = -0.554700196225;
   cblas_drotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, dbleps, "drotg(case 272)");
   gsl_test_rel(b, z_expected, dbleps, "drotg(case 273)");
   gsl_test_rel(c, c_expected, dbleps, "drotg(case 274)");
   gsl_test_rel(s, s_expected, dbleps, "drotg(case 275)");
  };


  {
   double a = -1.5;
   double b = 1.5;
   double c;
   double s;
   double r_expected = 2.12132034356;
   double z_expected = -1.41421356237;
   double c_expected = -0.707106781187;
   double s_expected = 0.707106781187;
   cblas_drotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, dbleps, "drotg(case 276)");
   gsl_test_rel(b, z_expected, dbleps, "drotg(case 277)");
   gsl_test_rel(c, c_expected, dbleps, "drotg(case 278)");
   gsl_test_rel(s, s_expected, dbleps, "drotg(case 279)");
  };


  {
   double a = -1;
   double b = -1.5;
   double c;
   double s;
   double r_expected = 1.80277563773;
   double z_expected = -1.80277563773;
   double c_expected = -0.554700196225;
   double s_expected = -0.832050294338;
   cblas_drotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, dbleps, "drotg(case 280)");
   gsl_test_rel(b, z_expected, dbleps, "drotg(case 281)");
   gsl_test_rel(c, c_expected, dbleps, "drotg(case 282)");
   gsl_test_rel(s, s_expected, dbleps, "drotg(case 283)");
  };


  {
   double a = -1;
   double b = -1;
   double c;
   double s;
   double r_expected = 1.41421356237;
   double z_expected = -1.41421356237;
   double c_expected = -0.707106781187;
   double s_expected = -0.707106781187;
   cblas_drotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, dbleps, "drotg(case 284)");
   gsl_test_rel(b, z_expected, dbleps, "drotg(case 285)");
   gsl_test_rel(c, c_expected, dbleps, "drotg(case 286)");
   gsl_test_rel(s, s_expected, dbleps, "drotg(case 287)");
  };


  {
   double a = -1;
   double b = -0.1;
   double c;
   double s;
   double r_expected = -1.00498756211;
   double z_expected = 0.099503719021;
   double c_expected = 0.99503719021;
   double s_expected = 0.099503719021;
   cblas_drotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, dbleps, "drotg(case 288)");
   gsl_test_rel(b, z_expected, dbleps, "drotg(case 289)");
   gsl_test_rel(c, c_expected, dbleps, "drotg(case 290)");
   gsl_test_rel(s, s_expected, dbleps, "drotg(case 291)");
  };


  {
   double a = -1;
   double b = 0;
   double c;
   double s;
   double r_expected = -1;
   double z_expected = 0;
   double c_expected = 1;
   double s_expected = 0;
   cblas_drotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, dbleps, "drotg(case 292)");
   gsl_test_rel(b, z_expected, dbleps, "drotg(case 293)");
   gsl_test_rel(c, c_expected, dbleps, "drotg(case 294)");
   gsl_test_rel(s, s_expected, dbleps, "drotg(case 295)");
  };


  {
   double a = -1;
   double b = 0.1;
   double c;
   double s;
   double r_expected = -1.00498756211;
   double z_expected = -0.099503719021;
   double c_expected = 0.99503719021;
   double s_expected = -0.099503719021;
   cblas_drotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, dbleps, "drotg(case 296)");
   gsl_test_rel(b, z_expected, dbleps, "drotg(case 297)");
   gsl_test_rel(c, c_expected, dbleps, "drotg(case 298)");
   gsl_test_rel(s, s_expected, dbleps, "drotg(case 299)");
  };


  {
   double a = -1;
   double b = 1;
   double c;
   double s;
   double r_expected = 1.41421356237;
   double z_expected = -1.41421356237;
   double c_expected = -0.707106781187;
   double s_expected = 0.707106781187;
   cblas_drotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, dbleps, "drotg(case 300)");
   gsl_test_rel(b, z_expected, dbleps, "drotg(case 301)");
   gsl_test_rel(c, c_expected, dbleps, "drotg(case 302)");
   gsl_test_rel(s, s_expected, dbleps, "drotg(case 303)");
  };


  {
   double a = -1;
   double b = 1.5;
   double c;
   double s;
   double r_expected = 1.80277563773;
   double z_expected = -1.80277563773;
   double c_expected = -0.554700196225;
   double s_expected = 0.832050294338;
   cblas_drotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, dbleps, "drotg(case 304)");
   gsl_test_rel(b, z_expected, dbleps, "drotg(case 305)");
   gsl_test_rel(c, c_expected, dbleps, "drotg(case 306)");
   gsl_test_rel(s, s_expected, dbleps, "drotg(case 307)");
  };


  {
   double a = -0.1;
   double b = -1.5;
   double c;
   double s;
   double r_expected = 1.50332963784;
   double z_expected = -15.0332963784;
   double c_expected = -0.0665190105238;
   double s_expected = -0.997785157857;
   cblas_drotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, dbleps, "drotg(case 308)");
   gsl_test_rel(b, z_expected, dbleps, "drotg(case 309)");
   gsl_test_rel(c, c_expected, dbleps, "drotg(case 310)");
   gsl_test_rel(s, s_expected, dbleps, "drotg(case 311)");
  };


  {
   double a = -0.1;
   double b = -1;
   double c;
   double s;
   double r_expected = 1.00498756211;
   double z_expected = -10.0498756211;
   double c_expected = -0.099503719021;
   double s_expected = -0.99503719021;
   cblas_drotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, dbleps, "drotg(case 312)");
   gsl_test_rel(b, z_expected, dbleps, "drotg(case 313)");
   gsl_test_rel(c, c_expected, dbleps, "drotg(case 314)");
   gsl_test_rel(s, s_expected, dbleps, "drotg(case 315)");
  };


  {
   double a = -0.1;
   double b = -0.1;
   double c;
   double s;
   double r_expected = 0.141421356237;
   double z_expected = -1.41421356237;
   double c_expected = -0.707106781187;
   double s_expected = -0.707106781187;
   cblas_drotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, dbleps, "drotg(case 316)");
   gsl_test_rel(b, z_expected, dbleps, "drotg(case 317)");
   gsl_test_rel(c, c_expected, dbleps, "drotg(case 318)");
   gsl_test_rel(s, s_expected, dbleps, "drotg(case 319)");
  };


  {
   double a = -0.1;
   double b = 0;
   double c;
   double s;
   double r_expected = -0.1;
   double z_expected = 0;
   double c_expected = 1;
   double s_expected = 0;
   cblas_drotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, dbleps, "drotg(case 320)");
   gsl_test_rel(b, z_expected, dbleps, "drotg(case 321)");
   gsl_test_rel(c, c_expected, dbleps, "drotg(case 322)");
   gsl_test_rel(s, s_expected, dbleps, "drotg(case 323)");
  };


  {
   double a = -0.1;
   double b = 0.1;
   double c;
   double s;
   double r_expected = 0.141421356237;
   double z_expected = -1.41421356237;
   double c_expected = -0.707106781187;
   double s_expected = 0.707106781187;
   cblas_drotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, dbleps, "drotg(case 324)");
   gsl_test_rel(b, z_expected, dbleps, "drotg(case 325)");
   gsl_test_rel(c, c_expected, dbleps, "drotg(case 326)");
   gsl_test_rel(s, s_expected, dbleps, "drotg(case 327)");
  };


  {
   double a = -0.1;
   double b = 1;
   double c;
   double s;
   double r_expected = 1.00498756211;
   double z_expected = -10.0498756211;
   double c_expected = -0.099503719021;
   double s_expected = 0.99503719021;
   cblas_drotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, dbleps, "drotg(case 328)");
   gsl_test_rel(b, z_expected, dbleps, "drotg(case 329)");
   gsl_test_rel(c, c_expected, dbleps, "drotg(case 330)");
   gsl_test_rel(s, s_expected, dbleps, "drotg(case 331)");
  };


  {
   double a = -0.1;
   double b = 1.5;
   double c;
   double s;
   double r_expected = 1.50332963784;
   double z_expected = -15.0332963784;
   double c_expected = -0.0665190105238;
   double s_expected = 0.997785157857;
   cblas_drotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, dbleps, "drotg(case 332)");
   gsl_test_rel(b, z_expected, dbleps, "drotg(case 333)");
   gsl_test_rel(c, c_expected, dbleps, "drotg(case 334)");
   gsl_test_rel(s, s_expected, dbleps, "drotg(case 335)");
  };


  {
   double a = 0;
   double b = -1.5;
   double c;
   double s;
   double r_expected = -1.5;
   double z_expected = 1;
   double c_expected = 0;
   double s_expected = 1;
   cblas_drotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, dbleps, "drotg(case 336)");
   gsl_test_rel(b, z_expected, dbleps, "drotg(case 337)");
   gsl_test_rel(c, c_expected, dbleps, "drotg(case 338)");
   gsl_test_rel(s, s_expected, dbleps, "drotg(case 339)");
  };


  {
   double a = 0;
   double b = -1;
   double c;
   double s;
   double r_expected = -1;
   double z_expected = 1;
   double c_expected = 0;
   double s_expected = 1;
   cblas_drotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, dbleps, "drotg(case 340)");
   gsl_test_rel(b, z_expected, dbleps, "drotg(case 341)");
   gsl_test_rel(c, c_expected, dbleps, "drotg(case 342)");
   gsl_test_rel(s, s_expected, dbleps, "drotg(case 343)");
  };


  {
   double a = 0;
   double b = -0.1;
   double c;
   double s;
   double r_expected = -0.1;
   double z_expected = 1;
   double c_expected = 0;
   double s_expected = 1;
   cblas_drotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, dbleps, "drotg(case 344)");
   gsl_test_rel(b, z_expected, dbleps, "drotg(case 345)");
   gsl_test_rel(c, c_expected, dbleps, "drotg(case 346)");
   gsl_test_rel(s, s_expected, dbleps, "drotg(case 347)");
  };


  {
   double a = 0;
   double b = 0;
   double c;
   double s;
   double r_expected = 0;
   double z_expected = 1;
   double c_expected = 1;
   double s_expected = 0;
   cblas_drotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, dbleps, "drotg(case 348)");
   gsl_test_rel(b, z_expected, dbleps, "drotg(case 349)");
   gsl_test_rel(c, c_expected, dbleps, "drotg(case 350)");
   gsl_test_rel(s, s_expected, dbleps, "drotg(case 351)");
  };


  {
   double a = 0;
   double b = 0.1;
   double c;
   double s;
   double r_expected = 0.1;
   double z_expected = 1;
   double c_expected = 0;
   double s_expected = 1;
   cblas_drotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, dbleps, "drotg(case 352)");
   gsl_test_rel(b, z_expected, dbleps, "drotg(case 353)");
   gsl_test_rel(c, c_expected, dbleps, "drotg(case 354)");
   gsl_test_rel(s, s_expected, dbleps, "drotg(case 355)");
  };


  {
   double a = 0;
   double b = 1;
   double c;
   double s;
   double r_expected = 1;
   double z_expected = 1;
   double c_expected = 0;
   double s_expected = 1;
   cblas_drotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, dbleps, "drotg(case 356)");
   gsl_test_rel(b, z_expected, dbleps, "drotg(case 357)");
   gsl_test_rel(c, c_expected, dbleps, "drotg(case 358)");
   gsl_test_rel(s, s_expected, dbleps, "drotg(case 359)");
  };


  {
   double a = 0;
   double b = 1.5;
   double c;
   double s;
   double r_expected = 1.5;
   double z_expected = 1;
   double c_expected = 0;
   double s_expected = 1;
   cblas_drotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, dbleps, "drotg(case 360)");
   gsl_test_rel(b, z_expected, dbleps, "drotg(case 361)");
   gsl_test_rel(c, c_expected, dbleps, "drotg(case 362)");
   gsl_test_rel(s, s_expected, dbleps, "drotg(case 363)");
  };


  {
   double a = 0.1;
   double b = -1.5;
   double c;
   double s;
   double r_expected = 1.50332963784;
   double z_expected = 15.0332963784;
   double c_expected = 0.0665190105238;
   double s_expected = -0.997785157857;
   cblas_drotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, dbleps, "drotg(case 364)");
   gsl_test_rel(b, z_expected, dbleps, "drotg(case 365)");
   gsl_test_rel(c, c_expected, dbleps, "drotg(case 366)");
   gsl_test_rel(s, s_expected, dbleps, "drotg(case 367)");
  };


  {
   double a = 0.1;
   double b = -1;
   double c;
   double s;
   double r_expected = 1.00498756211;
   double z_expected = 10.0498756211;
   double c_expected = 0.099503719021;
   double s_expected = -0.99503719021;
   cblas_drotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, dbleps, "drotg(case 368)");
   gsl_test_rel(b, z_expected, dbleps, "drotg(case 369)");
   gsl_test_rel(c, c_expected, dbleps, "drotg(case 370)");
   gsl_test_rel(s, s_expected, dbleps, "drotg(case 371)");
  };


  {
   double a = 0.1;
   double b = -0.1;
   double c;
   double s;
   double r_expected = 0.141421356237;
   double z_expected = 1.41421356237;
   double c_expected = 0.707106781187;
   double s_expected = -0.707106781187;
   cblas_drotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, dbleps, "drotg(case 372)");
   gsl_test_rel(b, z_expected, dbleps, "drotg(case 373)");
   gsl_test_rel(c, c_expected, dbleps, "drotg(case 374)");
   gsl_test_rel(s, s_expected, dbleps, "drotg(case 375)");
  };


  {
   double a = 0.1;
   double b = 0;
   double c;
   double s;
   double r_expected = 0.1;
   double z_expected = 0;
   double c_expected = 1;
   double s_expected = 0;
   cblas_drotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, dbleps, "drotg(case 376)");
   gsl_test_rel(b, z_expected, dbleps, "drotg(case 377)");
   gsl_test_rel(c, c_expected, dbleps, "drotg(case 378)");
   gsl_test_rel(s, s_expected, dbleps, "drotg(case 379)");
  };


  {
   double a = 0.1;
   double b = 0.1;
   double c;
   double s;
   double r_expected = 0.141421356237;
   double z_expected = 1.41421356237;
   double c_expected = 0.707106781187;
   double s_expected = 0.707106781187;
   cblas_drotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, dbleps, "drotg(case 380)");
   gsl_test_rel(b, z_expected, dbleps, "drotg(case 381)");
   gsl_test_rel(c, c_expected, dbleps, "drotg(case 382)");
   gsl_test_rel(s, s_expected, dbleps, "drotg(case 383)");
  };


  {
   double a = 0.1;
   double b = 1;
   double c;
   double s;
   double r_expected = 1.00498756211;
   double z_expected = 10.0498756211;
   double c_expected = 0.099503719021;
   double s_expected = 0.99503719021;
   cblas_drotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, dbleps, "drotg(case 384)");
   gsl_test_rel(b, z_expected, dbleps, "drotg(case 385)");
   gsl_test_rel(c, c_expected, dbleps, "drotg(case 386)");
   gsl_test_rel(s, s_expected, dbleps, "drotg(case 387)");
  };


  {
   double a = 0.1;
   double b = 1.5;
   double c;
   double s;
   double r_expected = 1.50332963784;
   double z_expected = 15.0332963784;
   double c_expected = 0.0665190105238;
   double s_expected = 0.997785157857;
   cblas_drotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, dbleps, "drotg(case 388)");
   gsl_test_rel(b, z_expected, dbleps, "drotg(case 389)");
   gsl_test_rel(c, c_expected, dbleps, "drotg(case 390)");
   gsl_test_rel(s, s_expected, dbleps, "drotg(case 391)");
  };


  {
   double a = 1;
   double b = -1.5;
   double c;
   double s;
   double r_expected = 1.80277563773;
   double z_expected = 1.80277563773;
   double c_expected = 0.554700196225;
   double s_expected = -0.832050294338;
   cblas_drotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, dbleps, "drotg(case 392)");
   gsl_test_rel(b, z_expected, dbleps, "drotg(case 393)");
   gsl_test_rel(c, c_expected, dbleps, "drotg(case 394)");
   gsl_test_rel(s, s_expected, dbleps, "drotg(case 395)");
  };


  {
   double a = 1;
   double b = -1;
   double c;
   double s;
   double r_expected = 1.41421356237;
   double z_expected = 1.41421356237;
   double c_expected = 0.707106781187;
   double s_expected = -0.707106781187;
   cblas_drotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, dbleps, "drotg(case 396)");
   gsl_test_rel(b, z_expected, dbleps, "drotg(case 397)");
   gsl_test_rel(c, c_expected, dbleps, "drotg(case 398)");
   gsl_test_rel(s, s_expected, dbleps, "drotg(case 399)");
  };


  {
   double a = 1;
   double b = -0.1;
   double c;
   double s;
   double r_expected = 1.00498756211;
   double z_expected = -0.099503719021;
   double c_expected = 0.99503719021;
   double s_expected = -0.099503719021;
   cblas_drotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, dbleps, "drotg(case 400)");
   gsl_test_rel(b, z_expected, dbleps, "drotg(case 401)");
   gsl_test_rel(c, c_expected, dbleps, "drotg(case 402)");
   gsl_test_rel(s, s_expected, dbleps, "drotg(case 403)");
  };


  {
   double a = 1;
   double b = 0;
   double c;
   double s;
   double r_expected = 1;
   double z_expected = 0;
   double c_expected = 1;
   double s_expected = 0;
   cblas_drotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, dbleps, "drotg(case 404)");
   gsl_test_rel(b, z_expected, dbleps, "drotg(case 405)");
   gsl_test_rel(c, c_expected, dbleps, "drotg(case 406)");
   gsl_test_rel(s, s_expected, dbleps, "drotg(case 407)");
  };


  {
   double a = 1;
   double b = 0.1;
   double c;
   double s;
   double r_expected = 1.00498756211;
   double z_expected = 0.099503719021;
   double c_expected = 0.99503719021;
   double s_expected = 0.099503719021;
   cblas_drotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, dbleps, "drotg(case 408)");
   gsl_test_rel(b, z_expected, dbleps, "drotg(case 409)");
   gsl_test_rel(c, c_expected, dbleps, "drotg(case 410)");
   gsl_test_rel(s, s_expected, dbleps, "drotg(case 411)");
  };


  {
   double a = 1;
   double b = 1;
   double c;
   double s;
   double r_expected = 1.41421356237;
   double z_expected = 1.41421356237;
   double c_expected = 0.707106781187;
   double s_expected = 0.707106781187;
   cblas_drotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, dbleps, "drotg(case 412)");
   gsl_test_rel(b, z_expected, dbleps, "drotg(case 413)");
   gsl_test_rel(c, c_expected, dbleps, "drotg(case 414)");
   gsl_test_rel(s, s_expected, dbleps, "drotg(case 415)");
  };


  {
   double a = 1;
   double b = 1.5;
   double c;
   double s;
   double r_expected = 1.80277563773;
   double z_expected = 1.80277563773;
   double c_expected = 0.554700196225;
   double s_expected = 0.832050294338;
   cblas_drotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, dbleps, "drotg(case 416)");
   gsl_test_rel(b, z_expected, dbleps, "drotg(case 417)");
   gsl_test_rel(c, c_expected, dbleps, "drotg(case 418)");
   gsl_test_rel(s, s_expected, dbleps, "drotg(case 419)");
  };


  {
   double a = 1.5;
   double b = -1.5;
   double c;
   double s;
   double r_expected = 2.12132034356;
   double z_expected = 1.41421356237;
   double c_expected = 0.707106781187;
   double s_expected = -0.707106781187;
   cblas_drotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, dbleps, "drotg(case 420)");
   gsl_test_rel(b, z_expected, dbleps, "drotg(case 421)");
   gsl_test_rel(c, c_expected, dbleps, "drotg(case 422)");
   gsl_test_rel(s, s_expected, dbleps, "drotg(case 423)");
  };


  {
   double a = 1.5;
   double b = -1;
   double c;
   double s;
   double r_expected = 1.80277563773;
   double z_expected = -0.554700196225;
   double c_expected = 0.832050294338;
   double s_expected = -0.554700196225;
   cblas_drotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, dbleps, "drotg(case 424)");
   gsl_test_rel(b, z_expected, dbleps, "drotg(case 425)");
   gsl_test_rel(c, c_expected, dbleps, "drotg(case 426)");
   gsl_test_rel(s, s_expected, dbleps, "drotg(case 427)");
  };


  {
   double a = 1.5;
   double b = -0.1;
   double c;
   double s;
   double r_expected = 1.50332963784;
   double z_expected = -0.0665190105238;
   double c_expected = 0.997785157857;
   double s_expected = -0.0665190105238;
   cblas_drotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, dbleps, "drotg(case 428)");
   gsl_test_rel(b, z_expected, dbleps, "drotg(case 429)");
   gsl_test_rel(c, c_expected, dbleps, "drotg(case 430)");
   gsl_test_rel(s, s_expected, dbleps, "drotg(case 431)");
  };


  {
   double a = 1.5;
   double b = 0;
   double c;
   double s;
   double r_expected = 1.5;
   double z_expected = 0;
   double c_expected = 1;
   double s_expected = 0;
   cblas_drotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, dbleps, "drotg(case 432)");
   gsl_test_rel(b, z_expected, dbleps, "drotg(case 433)");
   gsl_test_rel(c, c_expected, dbleps, "drotg(case 434)");
   gsl_test_rel(s, s_expected, dbleps, "drotg(case 435)");
  };


  {
   double a = 1.5;
   double b = 0.1;
   double c;
   double s;
   double r_expected = 1.50332963784;
   double z_expected = 0.0665190105238;
   double c_expected = 0.997785157857;
   double s_expected = 0.0665190105238;
   cblas_drotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, dbleps, "drotg(case 436)");
   gsl_test_rel(b, z_expected, dbleps, "drotg(case 437)");
   gsl_test_rel(c, c_expected, dbleps, "drotg(case 438)");
   gsl_test_rel(s, s_expected, dbleps, "drotg(case 439)");
  };


  {
   double a = 1.5;
   double b = 1;
   double c;
   double s;
   double r_expected = 1.80277563773;
   double z_expected = 0.554700196225;
   double c_expected = 0.832050294338;
   double s_expected = 0.554700196225;
   cblas_drotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, dbleps, "drotg(case 440)");
   gsl_test_rel(b, z_expected, dbleps, "drotg(case 441)");
   gsl_test_rel(c, c_expected, dbleps, "drotg(case 442)");
   gsl_test_rel(s, s_expected, dbleps, "drotg(case 443)");
  };


  {
   double a = 1.5;
   double b = 1.5;
   double c;
   double s;
   double r_expected = 2.12132034356;
   double z_expected = 1.41421356237;
   double c_expected = 0.707106781187;
   double s_expected = 0.707106781187;
   cblas_drotg(&a, &b, &c, &s);
   gsl_test_rel(a, r_expected, dbleps, "drotg(case 444)");
   gsl_test_rel(b, z_expected, dbleps, "drotg(case 445)");
   gsl_test_rel(c, c_expected, dbleps, "drotg(case 446)");
   gsl_test_rel(s, s_expected, dbleps, "drotg(case 447)");
  };


}
