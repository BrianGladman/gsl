/* Author: G. Jungman
 * RCS: $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_ellint.h"

#define gslMIN(x,y)  ((x) < (y) ? (x) : (y))
#define gslMAX(x,y)  ((x) > (y) ? (x) : (y))


/*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*/

static double inline gslMAX3(double x, double y, double z)
{
  double xy = gslMAX(x, y);
  return gslMAX(xy, z);
}
static double inline gslMAX4(double x, double y, double z, double w)
{
  double xy  = gslMAX(x,  y);
  double xyz = gslMAX(xy, z);
  return gslMAX(xyz, w);
}


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

/* based on Carlson's algorithms:
   [B. C. Carlson Numer. Math. 33 (1979)]
   
   see also:
   [B.C. Carlson, Special Functions of Applied Mathematics (1977)]
 */

/*
      DOUBLE PRECISION FUNCTION RC(X,Y,ERRTOL,IERR)
C
C          THIS FUNCTION SUBROUTINE COMPUTES THE ELEMENTARY INTEGRAL
C          RC(X,Y) = INTEGRAL FROM ZERO TO INFINITY OF
C
C                              -1/2     -1
C                    (1/2)(T+X)    (T+Y)  DT,
C
C          WHERE X IS NONNEGATIVE AND Y IS POSITIVE.  THE DUPLICATION
C          THEOREM IS ITERATED UNTIL THE VARIABLES ARE NEARLY EQUAL,
C          AND THE FUNCTION IS THEN EXPANDED IN TAYLOR SERIES TO FIFTH
C          ORDER.  LOGARITHMIC, INVERSE CIRCULAR, AND INVERSE HYPER-
C          BOLIC FUNCTIONS CAN BE EXPRESSED IN TERMS OF RC.  REFERENCE:
C          B. C. CARLSON, COMPUTING ELLIPTIC INTEGRALS BY DUPLICATION,
C          NUMER. MATH. 33 (1979), 1-16.  CODED BY B. C. CARLSON AND
C          ELAINE M. NOTIS, AMES LABORATORY-DOE, IOWA STATE UNIVERSITY,
C          AMES, IOWA 50011.  MARCH 1, 1980.
C
C          CHECK BY ADDITION THEOREM: RC(X,X+Z) + RC(Y,Y+Z) = RC(0,Z),
C          WHERE X, Y, AND Z ARE POSITIVE AND  X * Y = Z * Z.
C
C          LOLIM AND UPLIM DETERMINE THE RANGE OF VALID ARGUMENTS.
C          LOLIM IS NOT LESS THAN THE MACHINE MINIMUM MULTIPLIED BY 5.
C          UPLIM IS NOT GREATER THAN THE MACHINE MAXIMUM DIVIDED BY 5.
C
      DATA LOLIM/3.D-78/, UPLIM/1.D+75/
C
C          ON INPUT:
C
C          X AND Y ARE THE VARIABLES IN THE INTEGRAL RC(X,Y).
C
C          ERRTOL IS SET TO THE DESIRED ERROR TOLERANCE.
C          RELATIVE ERROR DUE TO TRUNCATION IS LESS THAN
C          16 * ERRTOL ** 6 / (1 - 2 * ERRTOL).
C
C          SAMPLE CHOICES:  ERRTOL   RELATIVE TRUNCATION
C                                    ERROR LESS THAN
C                           1.D-3    2.D-17
C                           3.D-3    2.D-14
C                           1.D-2    2.D-11
C                           3.D-2    2.D-8
C                           1.D-1    2.D-5
C
*/
int gsl_sf_ellint_RC_impl(double x, double y, double errtol, double * result)
{
  const double lolim = 5.0 * DBL_MIN + 1.0;
  const double uplim = 0.2 * DBL_MIN - 1.0;

  if(x < 0.0 || y < 0.0 || x + y < lolim) {
    return GSL_EDOM;
  }
  else if(gslMAX(x, y) < uplim) { 
    const double c1 = 1.0 / 7.0;
    const double c2 = 9.0 / 22.0;
    double xn = x;
    double yn = y;
    double mu, sn, lamda, s;
    while(1) {
      mu = (xn + yn + yn) / 3.0;
      sn = (yn + mu) / mu - 2.0;
      if (fabs(sn) < errtol) break;
      lamda = 2.0 * sqrt(xn) * sqrt(yn) + yn;
      xn = (xn + lamda) * 0.25;
      yn = (yn + lamda) * 0.25;
    }
    s = sn * sn * (0.3 + sn * (c1 + sn * (0.375 + sn * c2)));
    *result = (1.0 + s) / sqrt(mu);
    return GSL_SUCCESS;
  }
  else {
    return GSL_EDOM;
  }
}

/*
      DOUBLE PRECISION FUNCTION RD(X,Y,Z,ERRTOL,IERR)
C
C          THIS FUNCTION SUBROUTINE COMPUTES AN INCOMPLETE ELLIPTIC
C          INTEGRAL OF THE SECOND KIND,
C          RD(X,Y,Z) = INTEGRAL FROM ZERO TO INFINITY OF
C
C                                -1/2     -1/2     -3/2
C                      (3/2)(T+X)    (T+Y)    (T+Z)    DT,
C
C          WHERE X AND Y ARE NONNEGATIVE, X + Y IS POSITIVE, AND Z IS
C          POSITIVE.  IF X OR Y IS ZERO, THE INTEGRAL IS COMPLETE.
C          THE DUPLICATION THEOREM IS ITERATED UNTIL THE VARIABLES ARE
C          NEARLY EQUAL, AND THE FUNCTION IS THEN EXPANDED IN TAYLOR
C          SERIES TO FIFTH ORDER.  REFERENCE: B. C. CARLSON, COMPUTING
C          ELLIPTIC INTEGRALS BY DUPLICATION, NUMER. MATH. 33 (1979),
C          1-16.  CODED BY B. C. CARLSON AND ELAINE M. NOTIS, AMES
C          LABORATORY-DOE, IOWA STATE UNIVERSITY, AMES, IOWA 50011.
C          MARCH 1, 1980..
C
C          CHECK: RD(X,Y,Z) + RD(Y,Z,X) + RD(Z,X,Y)
C          = 3 / DSQRT(X * Y * Z), WHERE X, Y, AND Z ARE POSITIVE.
C
C
C          LOLIM AND UPLIM DETERMINE THE RANGE OF VALID ARGUMENTS.
C          LOLIM IS NOT LESS THAN 2 / (MACHINE MAXIMUM) ** (2/3).
C          UPLIM IS NOT GREATER THAN (0.1 * ERRTOL / MACHINE
C          MINIMUM) ** (2/3), WHERE ERRTOL IS DESCRIBED BELOW.
C          IN THE FOLLOWING TABLE IT IS ASSUMED THAT ERRTOL WILL
C          NEVER BE CHOSEN SMALLER THAN 1.D-5.
C
      DATA LOLIM/6.D-51/, UPLIM/1.D+48/
C
C          ERRTOL IS SET TO THE DESIRED ERROR TOLERANCE.
C          RELATIVE ERROR DUE TO TRUNCATION IS LESS THAN
C          3 * ERRTOL ** 6 / (1-ERRTOL) ** 3/2.
C
C          SAMPLE CHOICES:  ERRTOL   RELATIVE TRUNCATION
C                                    ERROR LESS THAN
C                           1.D-3    4.D-18
C                           3.D-3    3.D-15
C                           1.D-2    4.D-12
C                           3.D-2    3.D-9
C                           1.D-1    4.D-6
C
*/
int gsl_sf_ellint_RD_impl(double x, double y, double z, double errtol, double * result)
{
  const double lolim = 2.0/pow(DBL_MAX, 2./3.) + 1.0;
  const double uplim = pow(0.1*errtol/DBL_MIN, 2./3.) - 1.0;

  if(gslMIN(x,y) < 0.0 || gslMIN(x+y,z) < lolim) {
    return GSL_EDOM;
  }
  else if(gslMAX3(x,y,z) < uplim) {
    const double c1 = 3.0 / 14.0;
    const double c2 = 1.0 /  6.0;
    const double c3 = 9.0 / 22.0;
    const double c4 = 3.0 / 26.0;
    double xn = x;
    double yn = y;
    double zn = z;
    double sigma  = 0.0;
    double power4 = 1.0;
    double ea, eb, ec, ed, ef, s1, s2;
    double mu, xndev, yndev, zndev;
    while(1) {
      double xnroot, ynroot, znroot, lamda;
      double epslon;
      mu = (xn + yn + 3.0 * zn) * 0.2;
      xndev = (mu - xn) / mu;
      yndev = (mu - yn) / mu;
      zndev = (mu - zn) / mu;
      epslon = gslMAX3(fabs(xndev), fabs(yndev), fabs(zndev));
      if (epslon < errtol) break;
      xnroot = sqrt(xn);
      ynroot = sqrt(yn);
      znroot = sqrt(zn);
      lamda = xnroot * (ynroot + znroot) + ynroot * znroot;
      sigma  += power4 / (znroot * (zn + lamda));
      power4 *= 0.25;
      xn = (xn + lamda) * 0.25;
      yn = (yn + lamda) * 0.25;
      zn = (zn + lamda) * 0.25;
    }
    ea = xndev * yndev;
    eb = zndev * zndev;
    ec = ea - eb;
    ed = ea - 6.0 * eb;
    ef = ed + ec + ec;
    s1 = ed * (- c1 + 0.25 * c3 * ed - 1.5 * c4 * zndev * ef);
    s2 = zndev * (c2 * ef + zndev * (- c3 * ec + zndev * c4 * ea));
    *result = 3.0 * sigma + power4 * (1.0 + s1 + s2) / (mu * sqrt(mu));
    return GSL_SUCCESS;
  }
  else {
    return GSL_EDOM;
  }
}

/*
      DOUBLE PRECISION FUNCTION RF(X,Y,Z,ERRTOL,IERR)
C
C          THIS FUNCTION SUBROUTINE COMPUTES THE INCOMPLETE ELLIPTIC
C          INTEGRAL OF THE FIRST KIND,
C          RF(X,Y,Z) = INTEGRAL FROM ZERO TO INFINITY OF
C
C                                -1/2     -1/2     -1/2
C                      (1/2)(T+X)    (T+Y)    (T+Z)    DT,
C
C          WHERE X, Y, AND Z ARE NONNEGATIVE AND AT MOST ONE OF THEM
C          IS ZERO.  IF ONE OF THEM IS ZERO, THE INTEGRAL IS COMPLETE.
C          THE DUPLICATION THEOREM IS ITERATED UNTIL THE VARIABLES ARE
C          NEARLY EQUAL, AND THE FUNCTION IS THEN EXPANDED IN TAYLOR
C          SERIES TO FIFTH ORDER.  REFERENCE: B. C. CARLSON, COMPUTING
C          ELLIPTIC INTEGRALS BY DUPLICATION, NUMER. MATH. 33 (1979),
C          1-16.  CODED BY B. C. CARLSON AND ELAINE M. NOTIS, AMES
C          LABORATORY-DOE, IOWA STATE UNIVERSITY, AMES, IOWA 50011.
C          MARCH 1, 1980.
C
C          CHECK BY ADDITION THEOREM: RF(X,X+Z,X+W) + RF(Y,Y+Z,Y+W)
C          = RF(0,Z,W), WHERE X,Y,Z,W ARE POSITIVE AND X * Y = Z * W.
C
C          LOLIM AND UPLIM DETERMINE THE RANGE OF VALID ARGUMENTS.
C          LOLIM IS NOT LESS THAN THE MACHINE MINIMUM MULTIPLIED BY 5.
C          UPLIM IS NOT GREATER THAN THE MACHINE MAXIMUM DIVIDED BY 5.
C
      DATA LOLIM/3.D-78/, UPLIM/1.D+75/
C
C          ERRTOL IS SET TO THE DESIRED ERROR TOLERANCE.
C          RELATIVE ERROR DUE TO TRUNCATION IS LESS THAN
C          ERRTOL ** 6 / (4 * (1 - ERRTOL)).
C
C          SAMPLE CHOICES:  ERRTOL   RELATIVE TRUNCATION
C                                    ERROR LESS THAN
C                           1.D-3    3.D-19
C                           3.D-3    2.D-16
C                           1.D-2    3.D-13
C                           3.D-2    2.D-10
C                           1.D-1    3.D-7
*/
int gsl_sf_ellint_RF_impl(double x, double y, double z, double errtol, double * result)
{
  const double lolim = 5.0 * DBL_MIN + 1.0;
  const double uplim = 0.2 * DBL_MIN - 1.0;

  if(x < 0.0 || y < 0.0 || z < 0.0) {
    return GSL_EDOM;
  }
  else if(x+y < lolim || x + z < lolim || y+z < lolim) {
    return GSL_EDOM;
  }
  else if(gslMAX3(x,y,z) < uplim) { 
    const double c1 = 1.0 / 24.0;
    const double c2 = 3.0 / 44.0;
    const double c3 = 1.0 / 14.0;
    double xn = x;
    double yn = y;
    double zn = z;
    double mu, xndev, yndev, zndev, e2, e3, s;
    while(1) {
      double epslon, lamda;
      double xnroot, ynroot, znroot;
      mu = (xn + yn + zn) / 3.0;
      xndev = 2.0 - (mu + xn) / mu;
      yndev = 2.0 - (mu + yn) / mu;
      zndev = 2.0 - (mu + zn) / mu;
      epslon = gslMAX3(fabs(xndev), fabs(yndev), fabs(zndev));
      if (epslon < errtol) break;
      xnroot = sqrt(xn);
      ynroot = sqrt(yn);
      znroot = sqrt(zn);
      lamda = xnroot * (ynroot + znroot) + ynroot * znroot;
      xn = (xn + lamda) * 0.25;
      yn = (yn + lamda) * 0.25;
      zn = (zn + lamda) * 0.25;
    }
    e2 = xndev * yndev - zndev * zndev;
    e3 = xndev * yndev * zndev;
    s = 1.0 + (c1 * e2 - 0.10 - c2 * e3) * e2 + c3 * e3;
    *result = s / sqrt(mu);
    return GSL_SUCCESS;
  }
  else {
    return GSL_EDOM;
  }
}

/*
    DOUBLE PRECISION FUNCTION RJ(X,Y,Z,P,ERRTOL,IERR)
C
C          THIS FUNCTION SUBROUTINE COMPUTES AN INCOMPLETE ELLIPTIC
C          INTEGRAL OF THE THIRD KIND,
C          RJ(X,Y,Z,P) = INTEGRAL FROM ZERO TO INFINITY OF
C
C                                  -1/2     -1/2     -1/2     -1
C                        (3/2)(T+X)    (T+Y)    (T+Z)    (T+P)  DT,
C
C          WHERE X, Y, AND Z ARE NONNEGATIVE, AT MOST ONE OF THEM IS
C          ZERO, AND P IS POSITIVE.  IF X OR Y OR Z IS ZERO, THE
C          INTEGRAL IS COMPLETE.  THE DUPLICATION THEOREM IS ITERATED
C          UNTIL THE VARIABLES ARE NEARLY EQUAL, AND THE FUNCTION IS
C          THEN EXPANDED IN TAYLOR SERIES TO FIFTH ORDER.  REFERENCE:
C          B. C. CARLSON, COMPUTING ELLIPTIC INTEGRALS BY DUPLICATION,
C          NUMER. MATH. 33 (1979), 1-16.  CODED BY B. C. CARLSON AND
C          ELAINE M. NOTIS, AMES LABORATORY-DOE, IOWA STATE UNIVERSITY,
C          AMES, IOWA 50011.  MARCH 1, 1980.
C
C          CHECK BY ADDITION THEOREM: RJ(X,X+Z,X+W,X+P)
C          + RJ(Y,Y+Z,Y+W,Y+P) + (A-B) * RJ(A,B,B,A) + 3 / DSQRT(A)
C          = RJ(0,Z,W,P), WHERE X,Y,Z,W,P ARE POSITIVE AND X * Y
C          = Z * W,  A = P * P * (X+Y+Z+W),  B = P * (P+X) * (P+Y),
C          AND B - A = P * (P-Z) * (P-W).  THE SUM OF THE THIRD AND
C          FOURTH TERMS ON THE LEFT SIDE IS 3 * RC(A,B).
C
C          RC IS A FUNCTION COMPUTED BY AN EXTERNAL SUBROUTINE.
C
C          LOLIM AND UPLIM DETERMINE THE RANGE OF VALID ARGUMENTS.
C          LOLIM IS NOT LESS THAN THE CUBE ROOT OF THE VALUE
C          OF LOLIM USED IN THE SUBROUTINE FOR RC.
C          UPLIM IS NOT GREATER THAN 0.3 TIMES THE CUBE ROOT OF
C          THE VALUE OF UPLIM USED IN THE SUBROUTINE FOR RC.
C
      DATA LOLIM/2.D-26/, UPLIM/3.D+24/
C
C          ERRTOL IS SET TO THE DESIRED ERROR TOLERANCE.
C          RELATIVE ERROR DUE TO TRUNCATION OF THE SERIES FOR RJ
C          IS LESS THAN 3 * ERRTOL ** 6 / (1 - ERRTOL) ** 3/2.
C          AN ERROR TOLERANCE (ETOLRC) WILL BE PASSED TO THE SUBROUTINE
C          FOR RC TO MAKE THE TRUNCATION ERROR FOR RC LESS THAN FOR RJ.
C
C          SAMPLE CHOICES:  ERRTOL   RELATIVE TRUNCATION
C                                    ERROR LESS THAN
C                           1.D-3    4.D-18
C                           3.D-3    3.D-15
C                           1.D-2    4.D-12
C                           3.D-2    3.D-9
C                           1.D-1    4.D-6
C
*/
int gsl_sf_ellint_RJ_impl(double x, double y, double z, double p, double errtol, double * result)
{
  const double lolim =       pow(5.0 * DBL_MIN + 1.0, 1./3.);
  const double uplim = 0.3 * pow(0.2 * DBL_MIN - 1.0, 1./3.);

  if(x < 0.0 || y < 0.0 || y < 0.0) {
    return GSL_EDOM;
  }
  else if(x + y < lolim || x + z < lolim || y + z < lolim || p < lolim) {
    return GSL_EDOM;
  }
  else if(gslMAX4(x,y,z,p) < uplim) {
    const double c1 = 3.0 / 14.0;
    const double c2 = 1.0 /  3.0;
    const double c3 = 3.0 / 22.0;
    const double c4 = 3.0 / 26.0;
    double xn = x;
    double yn = y;
    double zn = z;
    double pn = p;
    double sigma = 0.0;
    double power4 = 1.0;
    double etolrc = 0.5 * errtol;
    double mu, xndev, yndev, zndev, pndev;
    double ea, eb, ec, e2, e3, s1, s2, s3;
    while(1) {
      double xnroot, ynroot, znroot;
      double lamda, alfa, beta;
      double epslon;
      double rcresult;
      int rcstatus;
      mu = (xn + yn + zn + pn + pn) * 0.2;
      xndev = (mu - xn) / mu;
      yndev = (mu - yn) / mu;
      zndev = (mu - zn) / mu;
      pndev = (mu - pn) / mu;
      epslon = gslmax4(fabs(xndev), fabs(yndev), fabs(zndev), fabs(pndev));
      if(epslon < errtol) break;
      xnroot = sqrt(xn);
      ynroot = sqrt(yn);
      znroot = sqrt(zn);
      lamda = xnroot * (ynroot + znroot) + ynroot * znroot;
      alfa = pn * (xnroot + ynroot + znroot) + xnroot * ynroot * znroot;
      alfa = alfa * alfa;
      beta = pn * (pn + lamda) * (pn + lamda);
      rcstatus = gsl_sf_ellint_RC_impl(alfa, beta, etolrc, &rcresult);
      if(rcstatus != GSL_SUCCESS) return GSL_EFAILED;
      sigma  += power4 * rcresult;
      power4 *= 0.25;
      xn = (xn + lamda) * 0.25;
      yn = (yn + lamda) * 0.25;
      zn = (zn + lamda) * 0.25;
      pn = (pn + lamda) * 0.25;
    }
    
    ea = xndev * (yndev + zndev) + yndev * zndev;
    eb = xndev * yndev * zndev;
    ec = pndev * pndev;
    e2 = ea - 3.0 * ec;
    e3 = eb + 2.0 * pndev * (ea - ec);
    s1 = 1.0 + e2 * (- c1 + 0.75 * c3 * e2 - 1.5 * c4 * e3);
    s2 = eb * (0.5 * c2 + pndev * (- c3 - c3 + pndev * c4));
    s3 = pndev * ea * (c2 - pndev * c3) - c2 * pndev * ec;
    *result = 3.0 * sigma + power4 * (s1 + s2 + s3) / (mu * sqrt(mu));
    return GSL_SUCCESS;
  }
  else {
    return GSL_EDOM;
  }
}

/* [Carlson, Numer. Math. 33 (1979) 1, (4.1)] */
int gsl_sf_ellint_F_impl(double phi, double k, double prec, double * result)
{
  double sin_phi  = sin(phi);
  double sin2_phi = sin_phi*sin_phi;
  double x = 1.0 - sin2_phi;
  double y = 1.0 - k*k*sin2_phi;
  double rf;
  int status = gsl_sf_ellint_RF_impl(x, y, 1.0, prec, &rf);
  if(status == GSL_SUCCESS) {
    *result = sin_phi * rf;
  }
  return status;
}

/* [Carlson, Numer. Math. 33 (1979) 1, (4.2)] */
int gsl_sf_ellint_E_impl(double phi, double k, double prec, double * result)
{
  double sin_phi  = sin(phi);
  double sin2_phi = sin_phi  * sin_phi;
  double sin3_phi = sin2_phi * sin_phi;
  double x = 1.0 - sin2_phi;
  double y = 1.0 - k*k*sin2_phi;
  double rf, rd;
  int rfstatus = gsl_sf_ellint_RF_impl(x, y, 1.0, prec, &rf);
  int rdstatus = gsl_sf_ellint_RD_impl(x, y, 1.0, prec, &rd);
  if(rfstatus == GSL_SUCCESS && rdstatus == GSL_SUCCESS) {
    *result = sin_phi * rf - k*k/3.0 * sin3_phi * rd;
    return GSL_SUCCESS;
  }
  else if(rfstatus == GSL_EDOM || rdstatus == GSL_EDOM) {
    return GSL_EDOM;
  }
  else {
    return GSL_EFAILED;
  }
}

/* [Carlson, Numer. Math. 33 (1979) 1, (4.3)] */
int gsl_sf_ellint_P_impl(double phi, double k, double n, double prec, double * result)
{
  double sin_phi  = sin(phi);
  double sin2_phi = sin_phi  * sin_phi;
  double sin3_phi = sin2_phi * sin_phi;
  double x = 1.0 - sin2_phi;
  double y = 1.0 - k*k*sin2_phi;
  double rf, rj;
  int rfstatus = gsl_sf_ellint_RF_impl(x, y, 1.0, prec, &rf);
  int rjstatus = gsl_sf_ellint_RJ_impl(x, y, 1.0, 1.0 + n*sin3_phi, prec, &rj);
  if(rfstatus == GSL_SUCCESS && rjstatus == GSL_SUCCESS) {
    *result = sin_phi * rf - n/3.0*k*k*sin3_phi * rj;
    return GSL_SUCCESS;
  }
  else if(rfstatus == GSL_EDOM || rjstatus == GSL_EDOM) {
    return GSL_EDOM;
  }
  else {
    return GSL_EFAILED;
  }
}

/* [Carlson, Numer. Math. 33 (1979) 1, (4.4)] */
int gsl_sf_ellint_D_impl(double phi, double k, double n, double prec, double * result)
{
  double sin_phi  = sin(phi);
  double sin2_phi = sin_phi  * sin_phi;
  double sin3_phi = sin2_phi * sin_phi;
  double x = 1.0 - sin2_phi;
  double y = 1.0 - k*k*sin2_phi;
  double rd;
  int status = gsl_sf_ellint_RD_impl(x, y, 1.0, prec, &rd);
  if(status == GSL_SUCCESS) {
    *result = sin3_phi/3.0 * rd;
  }
  return status;
}

/* [Carlson, Numer. Math. 33 (1979) 1, (4.5)] */
int gsl_sf_ellint_Kcomp_impl(double k, double prec, double *result)
{
  return gsl_sf_ellint_RF_impl(0.0, 1.0 - k*k, 1.0, prec, result);
}

/* [Carlson, Numer. Math. 33 (1979) 1, (4.6)] */
int gsl_sf_ellint_Ecomp_impl(double k, double prec, double * result)
{
  double y = 1.0 - k*k;
  double rf, rd;
  int rfstatus = gsl_sf_ellint_RF_impl(0.0, y, 1.0, prec, &rf);
  int rdstatus = gsl_sf_ellint_RF_impl(0.0, y, 1.0, prec, &rf);
  if(rfstatus == GSL_SUCCESS && rdstatus == GSL_SUCCESS) {
    *result = rf - k*k/3.0 * rd;
    return GSL_SUCCESS;
  }
  else if(rfstatus == GSL_EDOM || rdstatus == GSL_EDOM) {
    return GSL_EDOM;
  }
  else {
    return GSL_EFAILED;
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_ellint_Kcomp_e(double k, double prec, double * result)
{
  int status = gsl_sf_ellint_Kcomp_impl(k, prec, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_ellint_Kcomp_e", status);
  }
  return status;
}

int gsl_sf_ellint_Ecomp_e(double k, double prec, double * result)
{
  int status = gsl_sf_ellint_Ecomp_impl(k, prec, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_ellint_Ecomp_e", status);
  }
  return status;
}

int gsl_sf_ellint_F_e(double phi, double k, double prec, double * result)
{
  int status = gsl_sf_ellint_F_impl(phi, k, prec, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_ellint_F_e", status);
  }
  return status;
}

int gsl_sf_ellint_E_e(double phi, double k, double prec, double * result)
{
  int status = gsl_sf_ellint_E_impl(phi, k, prec, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_ellint_E_e", status);
  }
  return status;
}

int gsl_sf_ellint_P_e(double phi, double k, double n, double prec, double * result)
{
  int status = gsl_sf_ellint_P_impl(phi, k, n, prec, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_ellint_P_e", status);
  }
  return status;
}

int gsl_sf_ellint_D_e(double phi, double k, double n, double prec, double * result)
{
  int status = gsl_sf_ellint_D_impl(phi, k, n, prec, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_ellint_D_e", status);
  }
  return status;
}

int gsl_sf_ellint_RC_e(double x, double y, double prec, double * result)
{
  int status = gsl_sf_ellint_RC_impl(x, y, prec, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_ellint_RC_e", status);
  }
  return status;
}

int gsl_sf_ellint_RD_e(double x, double y, double z, double prec, double * result)
{
  int status = gsl_sf_ellint_RD_impl(x, y, z, prec, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_ellint_RD_e", status);
  }
  return status;
}

int gsl_sf_ellint_RF_e(double x, double y, double z, double prec, double * result)
{
  int status = gsl_sf_ellint_RF_impl(x, y, z, prec, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_ellint_RF_e", status);
  }
  return status;
}

int gsl_sf_ellint_RJ_e(double x, double y, double z, double p, double prec, double * result)
{
  int status = gsl_sf_ellint_RJ_impl(x, y, z, p, prec, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_ellint_RJ_e", status);
  }
  return status;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*-*/

double gsl_sf_ellint_Kcomp(double k, double prec)
{
  double y;
  int status = gsl_sf_ellint_Kcomp_impl(k, prec, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_ellint_Kcomp", status);
  }
  return y;
}

double gsl_sf_ellint_Ecomp(double k, double prec)
{
  double y;
  int status = gsl_sf_ellint_Ecomp_impl(k, prec, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_ellint_Ecomp", status);
  }
  return y;
}

double gsl_sf_ellint_F(double phi, double k, double prec)
{
  double y;
  int status = gsl_sf_ellint_F_impl(phi, k, prec, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_ellint_F", status);
  }
  return y;
}

double gsl_sf_ellint_E(double phi, double k, double prec)
{
  double y;
  int status = gsl_sf_ellint_E_impl(phi, k, prec, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_ellint_E", status);
  }
  return y;
}

double gsl_sf_ellint_P(double phi, double k, double n, double prec)
{
  double y;
  int status = gsl_sf_ellint_P_impl(phi, k, n, prec, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_ellint_P", status);
  }
  return y;
}

double gsl_sf_ellint_D(double phi, double k, double n, double prec)
{
  double y;
  int status = gsl_sf_ellint_D_impl(phi, k, n, prec, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_ellint_D", status);
  }
  return y;
}

double gsl_sf_ellint_RC(double x, double y, double prec)
{
  double yy;
  int status = gsl_sf_ellint_RC_impl(x, y, prec, &yy);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_ellint_RC", status);
  }
  return yy;
}

double gsl_sf_ellint_RD(double x, double y, double z, double prec)
{
  double yy;
  int status = gsl_sf_ellint_RD_impl(x, y, z, prec, &yy);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_ellint_RD", status);
  }
  return yy;
}

double gsl_sf_ellint_RF(double x, double y, double z, double prec)
{
  double yy;
  int status = gsl_sf_ellint_RF_impl(x, y, z, prec, &yy);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_ellint_RF", status);
  }
  return yy;
}

double gsl_sf_ellint_RJ(double x, double y, double z, double p, double prec)
{
  double yy;
  int status = gsl_sf_ellint_RJ_impl(x, y, z, p, prec, &yy);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_ellint_RJ", status);
  }
  return yy;
}
