#ifndef GSL_INTEGRATION_H
#define GSL_INTEGRATION_H

/* Quadrature of f() on (a,b) by Simpson rule.
   Attempts evaluation to precision specified by eps.
 */
double gsl_integ_simpson(double (*f)(double), double a, double b, double eps);

/* Table quadrature using the Simpson rule. */
double gsl_integ_simpson_table(const double * x, const double * y, int n);

/* Integrate a function of the form
 *
 *       f(x) / ((x-x0)^2 + w^2)
 *
 * where f(x) is assumed smoothly varying over the region
 * of integration. Uses change of variable and Simpson rule.
 */
double gsl_integ_lorenz(double (*f)(double),
			double x0, double w,
			double a, double b, double eps);


int gsl_integration_qk15 (double (*f) (double x),
			  double a, double b,
			  double * result, double * abserr,
			  double * resabs, double * resasc) ;

int gsl_integration_qk21 (double (*f) (double x),
			  double a, double b,
			  double * result, double * abserr,
			  double * resabs, double * resasc);

int gsl_integration_qk31 (double (*f) (double x),
			  double a, double b,
			  double * result, double * abserr,
			  double * resabs, double * resasc);

int gsl_integration_qk41 (double (*f) (double x),
			  double a, double b,
			  double * result, double * abserr,
			  double * resabs, double * resasc);

int gsl_integration_qk51 (double (*f) (double x),
			  double a, double b,
			  double * result, double * abserr,
			  double * resabs, double * resasc);

int gsl_integration_qk61 (double (*f) (double x),
			  double a, double b,
			  double * result, double * abserr,
			  double * resabs, double * resasc);

int gsl_integration_qk (const int n,
			const double xgk[], const double wg[], const double wgk[],
			double fv1[], double fv2[],
			double (*f) (double x),
			double a, double b,
			double * result, double * abserr,
			double * resabs, double * resasc) ;

int gsl_integration_qng (double (*f) (double x),
			 double a, double b,
			 double epsabs, double epsrel,
			 double * result, double * abserr,
			 int * neval);


#endif /* GSL_INTEGRATION_H */
