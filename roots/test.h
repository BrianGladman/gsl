/* $Id$ */

void
test_macros();

int
test_bisection (const char * description,
		double (* f)(double), 
		double lower_bound, double upper_bound, 
		double correct_root);

int
test_bisection_failure (const char * description,
			double (* f)(double), 
			double lower_bound, double upper_bound, 
			double correct_root);

int
test_falsepos (const char * description,
	       double (* f)(double), 
	       double lower_bound, double upper_bound, 
	       double correct_root);

int
test_falsepos_failure (const char * description,
		       double (* f)(double), 
		       double lower_bound, double upper_bound, 
		       double correct_root);
int
test_secant (const char * description,
	     double (* f)(double), 
	     double guess1, double guess2,
	     double correct_root);

int
test_secant_failure (const char * description,
		     double (* f)(double), 
		     double guess1, double guess2,
		     double correct_root);

int
test_newton(const char * description,
	    void (* fdf)(double *, double *, double, int, int),
	    double guess, double correct_root);

int
test_newton_failure(const char * description,
		    void (* fdf)(double *, double *, double, int, int),
		    double guess, double correct_root);

void
usage(void);

void
error_handler(const char * reason, const char * file, int line);

double
test_hairy_1(double x);

void
test_hairy_1_fdf(double * y, double * yprime, double x, int y_wanted,
                          int yprime_wanted);

double
test_hairy_2(double x);

void
test_hairy_2_fdf(double * y, double * yprime, double x, int y_wanted,
                          int yprime_wanted);

double
test_hairy_3(double x);

void
test_hairy_3_fdf(double * y, double * yprime, double x, int y_wanted,
                          int yprime_wanted);

double
test_hairy_4(double x);

void
test_hairy_4_fdf(double * y, double * yprime, double x, int y_wanted,
                          int yprime_wanted);

double
test_hairy_5(double x);

void
test_hairy_5_fdf(double * y, double * yprime, double x, int y_wanted,
                          int yprime_wanted);

double
test_hairy_6(double x);

void
test_hairy_6_fdf(double * y, double * yprime, double x, int y_wanted,
                          int yprime_wanted);

void
sin_fdf(double * y, double * yprime, double x, int y_wanted,
        int yprime_wanted);

void
cos_fdf(double * y, double * yprime, double x, int y_wanted,
        int yprime_wanted);

