
/* $Id$ */

void
  test_macros (void);

void
  test_bisection (const char *description,
		  double (*f) (double),
		  double lower_bound, double upper_bound,
		  double correct_root);

void
  test_bisection_failure (const char *description,
			  double (*f) (double),
			  double lower_bound, double upper_bound,
			  double correct_root);

void
  test_brent (const char *description,
		  double (*f) (double),
		  double lower_bound, double upper_bound,
		  double correct_root);

void
  test_brent_failure (const char *description,
			  double (*f) (double),
			  double lower_bound, double upper_bound,
			  double correct_root);

void
  test_falsepos (const char *description,
		 double (*f) (double),
		 double lower_bound, double upper_bound,
		 double correct_root);

void
  test_falsepos_failure (const char *description,
			 double (*f) (double),
			 double lower_bound, double upper_bound,
			 double correct_root);
void
  test_secant (const char *description,
	       double (*f) (double),
	       double guess1, double guess2,
	       double correct_root);

void
  test_secant_failure (const char *description,
		       double (*f) (double),
		       double guess1, double guess2,
		       double correct_root);

void
  test_newton (const char *description,
	       double (*f) (double),
	       double (*df) (double),
	       void (*fdf) (double, double *, double *),
	       double guess, double correct_root);

void
  test_newton_failure (const char *description,
		       double (*f) (double),
		       double (*df) (double),
		       void (*fdf) (double, double *, double *),
		       double guess, double correct_root);

void
  usage (void);

void
  error_handler (const char *reason, const char *file, int line);

double
  func1 (double x);

double
  func1_df (double x);

void
  func1_fdf (double x, double *y, double *yprime);

double
  func2 (double x);

double
  func2_df (double x);

void
  func2_fdf (double x, double *y, double *yprime);

double
  func3 (double x);

double
  func3_df (double x);

void
  func3_fdf (double x, double *y, double *yprime);

double
  func4 (double x);

double
  func4_df (double x);

void
  func4_fdf (double x, double *y, double *yprime);

double
  func5 (double x);

double
  func5_df (double x);

void
  func5_fdf (double x, double *y, double *yprime);

double
  func6 (double x);

double
  func6_df (double x);

void
  func6_fdf (double x, double *y, double *yprime);

double
  sin_f (double x);

double
  sin_df (double x);

void
  sin_fdf (double x, double *y, double *yprime);

double
  cos_f (double x);

double
  cos_df (double x);

void
  cos_fdf (double x, double *y, double *yprime);
