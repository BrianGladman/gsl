void rosenbrock_initpt (gsl_vector * x);
int rosenbrock_f (const gsl_vector * x, void *params, gsl_vector * f);
int rosenbrock_df (const gsl_vector * x, void *params, gsl_matrix * df);
int rosenbrock_fdf (const gsl_vector * x, void *params,	gsl_vector * f, gsl_matrix * df);

