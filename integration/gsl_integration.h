#ifndef GSL_INTEGRATION_H
#define GSL_INTEGRATION_H
#include <stdlib.h>

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


typedef void gsl_integration_rule_t (double (*f)(double x), 
				     const double a, const double b,
				     double * result, double * abserr,
				     double * defabs, double * resabs) ;
       
void gsl_integration_qk15 (double (*f) (double x),
			  const double a, const double b,
			  double * result, double * abserr,
			  double * resabs, double * resasc) ;

void gsl_integration_qk21 (double (*f) (double x),
			  const double a, const double b,
			  double * result, double * abserr,
			  double * resabs, double * resasc);

void gsl_integration_qk31 (double (*f) (double x),
			  double a, double b,
			  double * result, double * abserr,
			  double * resabs, double * resasc);

void gsl_integration_qk41 (double (*f) (double x),
			  double a, double b,
			  double * result, double * abserr,
			  double * resabs, double * resasc);

void gsl_integration_qk51 (double (*f) (double x),
			  double a, double b,
			  double * result, double * abserr,
			  double * resabs, double * resasc);

void gsl_integration_qk61 (double (*f) (double x),
			  const double a, const double b,
			  double * result, double * abserr,
			  double * resabs, double * resasc);

void gsl_integration_qk (const int n,
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
			 size_t * neval);

int
gsl_integration_qage (double (*f)(double x),
		      double a, double b,
		      double epsabs, double epsrel,
		      int key,
		      size_t limit,
		      double alist[], double blist[], double rlist[], 
		      double elist[], size_t iord[], size_t * last,
		      double * result, double * abserr, size_t * neval) ;

int
gsl_integration_qage_impl (double (*f)(double x),
			   const double a, const double b,
			   const double epsabs, const double epsrel,
			   const size_t limit,
			   double alist[], double blist[], double rlist[], 
			   double elist[], size_t iord[], size_t * last,
			   double * result, double * abserr, size_t * nqeval,
			   gsl_integration_rule_t * const q) ;



/* The low-level integration rules in QUADPACK are identified by small
   integers (1-6). We'll use symbolic constants to refer to them. 

   Don't change the values 1-6, we need those to compute the number of
   function evaluations used by the rule */

enum {   
  GSL_INTEG_GAUSS15 = 1,  /* 15 point Gauss-Kronrod rule */
  GSL_INTEG_GAUSS21 = 2,  /* 21 point Gauss-Kronrod rule */
  GSL_INTEG_GAUSS31 = 3,  /* 31 point Gauss-Kronrod rule */
  GSL_INTEG_GAUSS41 = 4,  /* 41 point Gauss-Kronrod rule */
  GSL_INTEG_GAUSS51 = 5,  /* 51 point Gauss-Kronrod rule */
  GSL_INTEG_GAUSS61 = 6   /* 61 point Gauss-Kronrod rule */
} ;

#endif /* GSL_INTEGRATION_H */
