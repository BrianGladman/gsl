void rosenbrock_initpt (gsl_vector * x);
int rosenbrock_f (const gsl_vector * x, void *params, gsl_vector * f);
int rosenbrock_df (const gsl_vector * x, void *params, gsl_matrix * df);
int rosenbrock_fdf (const gsl_vector * x, void *params,	gsl_vector * f, gsl_matrix * df);

void powellsingular_initpt (gsl_vector * x);
int powellsingular_f (const gsl_vector * x, void *params, gsl_vector * f);
int powellsingular_df (const gsl_vector * x, void *params, gsl_matrix * df);
int powellsingular_fdf (const gsl_vector * x, void *params,	gsl_vector * f, gsl_matrix * df);

