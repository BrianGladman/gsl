/* $Id$ */
#ifndef _GSL_SPECIAL_H_
#define _GSL_SPECIAL_H_

double GSL_erf(double x);
double GSL_erfc(double x);

double GSL_erfc_asymptotic(double x);
double GSL_erfseries(double x);

double GSL_log_erfc(double x);
double GSL_log_erfc_asymptotic(double x);

double GSL_Q(double x);
double GSL_Z(double x);

#endif /* _GSL_SPECIAL_H_ */
