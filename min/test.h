typedef double simple_function (double x);
typedef struct { simple_function * f; simple_function * df; } function_pair ;

gsl_function create_function (simple_function * f) ;
double eval_function (double x, void * params) ;

void
test_f_e (const gsl_min_fminimizer_type * T, 
	  const char * description, gsl_function *f,
	  double lower_bound, double minimum, double upper_bound, 
          double correct_minimum);

void
test_f (const gsl_min_fminimizer_type * T, 
        const char * description, gsl_function *f,
	double lower_bound, double middle, double upper_bound, 
        double correct_minimum);

int
test_bracket (const char * description,gsl_function *f,double lower_bound, 
	      double upper_bound, int max);

double func1 (double x);
double func2 (double x);
double func3 (double x);
