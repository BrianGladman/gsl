
typedef void (*initpt_function) (gsl_vector * x);

extern gsl_multiroot_function_fdf rosenbrock;
void rosenbrock_initpt (gsl_vector * x);
int rosenbrock_f (const gsl_vector * x, void *params, gsl_vector * f);
int rosenbrock_df (const gsl_vector * x, void *params, gsl_matrix * df);
int rosenbrock_fdf (const gsl_vector * x, void *params,	gsl_vector * f, gsl_matrix * df);

extern gsl_multiroot_function_fdf powellsing;
void powellsing_initpt (gsl_vector * x);
int powellsing_f (const gsl_vector * x, void *params, gsl_vector * f);
int powellsing_df (const gsl_vector * x, void *params, gsl_matrix * df);
int powellsing_fdf (const gsl_vector * x, void *params, gsl_vector * f, gsl_matrix * df);

extern gsl_multiroot_function_fdf powellscal;
void powellscal_initpt (gsl_vector * x);
int powellscal_f (const gsl_vector * x, void *params, gsl_vector * f);
int powellscal_df (const gsl_vector * x, void *params, gsl_matrix * df);
int powellscal_fdf (const gsl_vector * x, void *params, gsl_vector * f, gsl_matrix * df);

extern gsl_multiroot_function_fdf wood;
void wood_initpt (gsl_vector * x);
int wood_f (const gsl_vector * x, void *params, gsl_vector * f);
int wood_df (const gsl_vector * x, void *params, gsl_matrix * df);
int wood_fdf (const gsl_vector * x, void *params, gsl_vector * f, gsl_matrix * df);

extern gsl_multiroot_function_fdf helical;
void helical_initpt (gsl_vector * x);
int helical_f (const gsl_vector * x, void *params, gsl_vector * f);
int helical_df (const gsl_vector * x, void *params, gsl_matrix * df);
int helical_fdf (const gsl_vector * x, void *params, gsl_vector * f, gsl_matrix * df);
