#ifndef __ROOTS_H__
#define __ROOTS_H__

int
gsl_root_bisection (double * root, double (* f)(double), double * lower_bound,
                    double * upper_bound, double epsilon,
                    unsigned int max_iterations);

#endif /* __ROOTS_H__ */
