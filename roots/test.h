/* $Id$ */

void
gsl_test_root_macros(int * num_passed, int * num_failed);

void
gsl_test_root_bisection(double (* f)(double), double lower_bound,
                        double upper_bound, double cor_root, int * num_passed,
                        int * num_failed);

void
gsl_test_root_falsepos(double (* f)(double), double lower_bound,
                       double upper_bound, double cor_root, int * num_passed,
                       int * num_failed);

void
error_handler(const char * reason, const char * file, int line);

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
