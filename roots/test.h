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
gsl_test_root_secant(double (* f)(double), double guess1, double guess2,
                     double cor_root, int * num_passed, int * num_failed,
                     int failure_desired);

void
gsl_test_root_newton(void (* fdf)(double *, double *, double, int, int),
                     double guess, double cor_root, int * num_passed,
                     int * num_failed, int failure_desired);

void
print_usage(void);

void
error_handler(const char * reason, const char * file, int line);

double
gsl_test_root_hairy_1(double x);

void
gsl_test_root_hairy_1_fdf(double * y, double * yprime, double x, int y_wanted,
                          int yprime_wanted);

double
gsl_test_root_hairy_2(double x);

void
gsl_test_root_hairy_2_fdf(double * y, double * yprime, double x, int y_wanted,
                          int yprime_wanted);

double
gsl_test_root_hairy_3(double x);

void
gsl_test_root_hairy_3_fdf(double * y, double * yprime, double x, int y_wanted,
                          int yprime_wanted);

double
gsl_test_root_hairy_4(double x);

void
gsl_test_root_hairy_4_fdf(double * y, double * yprime, double x, int y_wanted,
                          int yprime_wanted);

double
gsl_test_root_hairy_5(double x);

void
gsl_test_root_hairy_5_fdf(double * y, double * yprime, double x, int y_wanted,
                          int yprime_wanted);

double
gsl_test_root_hairy_6(double x);

void
gsl_test_root_hairy_6_fdf(double * y, double * yprime, double x, int y_wanted,
                          int yprime_wanted);

void
sin_fdf(double * y, double * yprime, double x, int y_wanted,
        int yprime_wanted);

void
cos_fdf(double * y, double * yprime, double x, int y_wanted,
        int yprime_wanted);

