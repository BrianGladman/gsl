
/* $Id$ */

typedef double simple_function (double x);
typedef struct { simple_function * f; simple_function * df; } function_pair ;

gsl_function create_function (simple_function * f) ;
double eval_function (double x, void * params) ;

gsl_fdf create_fdf (simple_function * f, simple_function * df);
double eval_fdf_f (double x, void * params);
double eval_fdf_df (double x, void * params);
void eval_fdf (double x, void * params, double * y1, double * y2);

void
  test_macros (void);

void
  test_roots (void);

void
  test_poly (void);

void
test_f (const gsl_root_f_solver_type * T, 
	const char * description, gsl_function *f,
	double lower_bound, double upper_bound, double correct_root);

void
test_f_e (const gsl_root_f_solver_type * T, const char * description, 
	  gsl_function *f,
	  double lower_bound, double upper_bound, double correct_root);

void
test_fdf (const gsl_root_fdf_solver_type * T, const char * description, 
	  gsl_fdf *fdf, double root, double correct_root);

void
test_fdf_e (const gsl_root_fdf_solver_type * T, const char * description, 
	    gsl_fdf *fdf, double root, double correct_root);


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
