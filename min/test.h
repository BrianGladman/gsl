typedef double simple_function (double x);
typedef struct { simple_function * f; simple_function * df; } function_pair ;

gsl_function create_function (simple_function * f) ;
double eval_function (double x, void * params) ;

void
test_f (const gsl_min_fsolver_type * T, 
        const char * description, gsl_function *f,
	double lower_bound, double middle, double upper_bound, 
        double correct_minimum);
