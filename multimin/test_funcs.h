
typedef void (*initpt_function) (gsl_vector * x);

extern gsl_multimin_function_fdf rosenbrock;
void rosenbrock_initpt (gsl_vector * x);
double rosenbrock_f (const gsl_vector * x, void *params);
void rosenbrock_df (const gsl_vector * x, void *params, gsl_vector * df);
void rosenbrock_fdf (const gsl_vector * x, void *params, double * f, gsl_vector * df);

extern gsl_multimin_function_fdf wood;
void wood_initpt (gsl_vector * x);
double wood_f (const gsl_vector * x, void *params);
void wood_df (const gsl_vector * x, void *params, gsl_vector * df);
void wood_fdf (const gsl_vector * x, void *params, double * f, gsl_vector * df);

extern gsl_multimin_function_fdf roth;
void roth_initpt (gsl_vector * x);
double roth_f (const gsl_vector * x, void *params);
void roth_df (const gsl_vector * x, void *params, gsl_vector * df);
void roth_fdf (const gsl_vector * x, void *params, double * f, gsl_vector * df);

