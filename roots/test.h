/* $Id$ */

void
gsl_test_root_macros(int * num_passed, int * num_failed);

void
gsl_test_root_bisection(double (* f)(double), double lower_bound,
                        double upper_bound, double cur_root, int * num_passed,
                        int * num_failed);

double
gsl_test_root_hairy_1(double x);

double
gsl_test_root_hairy_2(double x);

double
gsl_test_root_hairy_3(double x);

double
gsl_test_root_hairy_4(double x);

double
gsl_test_root_hairy_5(double x);
