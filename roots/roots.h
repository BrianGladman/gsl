#ifndef __ROOTS_H__
#define __ROOTS_H__
double gsl_newton1D(double (*fn)(double x), double (*dfn)(double x),
		    double guess);
#endif /* __ROOTS_H__ */
