#ifndef GSL_INTEGRATION_H
#define GSL_INTEGRATION_H
#include <stdlib.h>
#include <gsl_math.h>

typedef struct {
  size_t limit ;
  double * alist ;
  double * blist ;
  double * rlist ;
  double * elist ;
  size_t * iord ;
} gsl_integration_workspace ;

typedef struct {
  size_t npts ;
  unsigned int * level;
  unsigned int * ndin ;
} gsl_integration_workspace_pts ;

gsl_integration_workspace * 
gsl_integration_workspace_alloc (size_t n) ;

void
gsl_integration_workspace_free (gsl_integration_workspace * w) ;

gsl_integration_workspace_pts *
gsl_integration_workspace_pts_alloc (const size_t npts);

void
gsl_integration_workspace_pts_free (gsl_integration_workspace_pts * w);

typedef void gsl_integration_rule_t (const gsl_function *f, 
				     const double a, const double b,
				     double * result, double * abserr,
				     double * defabs, double * resabs) ;
       
void gsl_integration_qk15 (const gsl_function *f,
			  const double a, const double b,
			  double * result, double * abserr,
			  double * resabs, double * resasc) ;

void gsl_integration_qk21 (const gsl_function *f,
			  const double a, const double b,
			  double * result, double * abserr,
			  double * resabs, double * resasc);

void gsl_integration_qk31 (const gsl_function *f,
			  double a, double b,
			  double * result, double * abserr,
			  double * resabs, double * resasc);

void gsl_integration_qk41 (const gsl_function *f,
			  double a, double b,
			  double * result, double * abserr,
			  double * resabs, double * resasc);

void gsl_integration_qk51 (const gsl_function *f,
			  double a, double b,
			  double * result, double * abserr,
			  double * resabs, double * resasc);

void gsl_integration_qk61 (const gsl_function *f,
			  const double a, const double b,
			  double * result, double * abserr,
			  double * resabs, double * resasc);

void gsl_integration_qk (const int n,
			 const double xgk[], const double wg[], const double wgk[],
			 double fv1[], double fv2[],
			 const gsl_function *f,
			 double a, double b,
			 double * result, double * abserr,
			 double * resabs, double * resasc) ;

int gsl_integration_qng (const gsl_function *f,
			 double a, double b,
			 double epsabs, double epsrel,
			 double * result, double * abserr,
			 size_t * neval);

int
gsl_integration_qag (const gsl_function *f,
		      double a, double b,
		      double epsabs, double epsrel,
		      int key,
		      gsl_integration_workspace * workspace,
		      size_t * last,
		      double * result, double * abserr, size_t * neval) ;

int
gsl_integration_qag_impl (const gsl_function *f,
			   const double a, const double b,
			   const double epsabs, const double epsrel,
			   gsl_integration_workspace * workspace,
			   size_t * last,
			   double * result, double * abserr, size_t * nqeval,
			   gsl_integration_rule_t * const q) ;

int
gsl_integration_qagi (gsl_function *f,
		      double epsabs, double epsrel,
		      gsl_integration_workspace * workspace,
		      size_t * last,
		      double * result, double * abserr, size_t * neval) ;

int
gsl_integration_qagiu (gsl_function *f,
		       double a,
		       double epsabs, double epsrel,
		       gsl_integration_workspace * workspace,
		       size_t * last,
		       double * result, double * abserr, size_t * neval) ;

int
gsl_integration_qagil (gsl_function *f,
		       double b,
		       double epsabs, double epsrel,
		       gsl_integration_workspace * workspace,
		       size_t * last,
		       double * result, double * abserr, size_t * neval) ;


void gsl_integration_qki (const int n,
			 const double xgk[], const double wg[], const double wgk[],
			 double fv1[], double fv2[],
			 const gsl_function *f,
			 double a, double b,
			 double * result, double * abserr,
			 double * resabs, double * resasc) ;

void gsl_integration_qk15i (const gsl_function *f,
			  const double a, const double b,
			  double * result, double * abserr,
			  double * resabs, double * resasc) ;


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

void
gsl_integration_qelg (size_t * n, double epstab[], 
		      double * result, double * abserr,
		      double res3la[], size_t * nres) ;

int
gsl_integration_qags (const gsl_function *f,
		      double a, double b,
		      double epsabs, double epsrel,
		      gsl_integration_workspace * workspace, size_t * last,
		      double * result, double * abserr, size_t * neval) ;

int
gsl_integration_qags_impl (const gsl_function *f, 
			   double a, double b, 
			   double epsabs, double epsrel,
			   gsl_integration_workspace * workspace,
			   double * result, double * abserr, 
			   size_t * last, size_t * nqeval,
			   gsl_integration_rule_t * q) ;
int
gsl_integration_qagp (const gsl_function *f,
		      double * pts, size_t npts,
		      double epsabs, double epsrel,
		      gsl_integration_workspace * workspace,
		      gsl_integration_workspace_pts * workspace_pts,
		      size_t * last,
		      double * result, double * abserr, size_t * neval);

int
gsl_integration_qagp_impl (const gsl_function *f,
			   const double *pts, const size_t npts,
			   double epsabs, double epsrel,
			   gsl_integration_workspace * workspace,
			   gsl_integration_workspace_pts * workspace_pts,
			   double *result, double *abserr,
			   size_t * last, size_t * nqeval,
			   gsl_integration_rule_t * const q);

#endif /* GSL_INTEGRATION_H */
