#ifndef __ROOTS_H__
#define __ROOTS_H__
void gsl_set_newton_epsilon(double new_val);
double gsl_get_newton_epsilon();
double gsl_newton1D(double (*fn)(double x), double (*dfn)(double x),
		    double guess);
#endif /* __ROOTS_H__ */
