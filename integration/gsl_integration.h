#ifndef GSL_INTEGRATION_H_
#define GSL_INTEGRATION_H_

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


#endif /* !GSL_INTEGRATION_H_ */
