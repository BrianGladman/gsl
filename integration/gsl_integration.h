#ifndef GSL_INTEGRATION_H
#define GSL_INTEGRATION_H
#include <stdlib.h>
#include <gsl/gsl_math.h>

/* Workspace for adaptive integrators */

typedef struct
  {
    size_t limit;
    size_t size;
    size_t nrmax;
    size_t i;
    size_t maximum_level;
    double *alist;
    double *blist;
    double *rlist;
    double *elist;
    size_t *order;
    size_t *level;
  }
gsl_integration_workspace;

gsl_integration_workspace *
  gsl_integration_workspace_alloc (size_t n);

void
  gsl_integration_workspace_free (gsl_integration_workspace * w);


/* Workspace for QAWS integrator */

typedef struct
{
  double alpha;
  double beta;
  int mu;
  int nu;
  double ri[25];
  double rj[25];
  double rg[25];
  double rh[25];
}
gsl_integration_qaws_table;

gsl_integration_qaws_table * 
gsl_integration_qaws_table_alloc (double alpha, double beta, int mu, int nu);

int
gsl_integration_qaws_table_set (gsl_integration_qaws_table * t,
                                double alpha, double beta, int mu, int nu);

void
gsl_integration_qaws_table_free (gsl_integration_qaws_table * t);

/* Workspace for QAWO integrator */

enum gsl_integration_qawo_enum { GSL_INTEG_COSINE, GSL_INTEG_SINE };

typedef struct
{
  size_t i;
  size_t n;
  double omega;
  double L;
  double par;
  enum gsl_integration_qawo_enum sine;
  double *chebmo;
}
gsl_integration_qawo_table;

gsl_integration_qawo_table * 
gsl_integration_qawo_table_alloc (double omega, double L, 
				      enum gsl_integration_qawo_enum sine,
				      size_t n);

int
gsl_integration_qawo_table_set (gsl_integration_qawo_table * t,
				    double omega, double L,
				    enum gsl_integration_qawo_enum sine);

int
gsl_integration_qawo_table_set_length (gsl_integration_qawo_table * t,
					   double L);

void
gsl_integration_qawo_table_free (gsl_integration_qawo_table * t);


/* Definition of an integration rule */

typedef void gsl_integration_rule (const gsl_function * f,
				   double a, double b,
				   double *result, double *abserr,
				   double *defabs, double *resabs);

void gsl_integration_qk15 (const gsl_function * f, double a, double b,
			   double *result, double *abserr,
			   double *resabs, double *resasc);

void gsl_integration_qk21 (const gsl_function * f, double a, double b,
			   double *result, double *abserr,
			   double *resabs, double *resasc);

void gsl_integration_qk31 (const gsl_function * f, double a, double b,
			   double *result, double *abserr,
			   double *resabs, double *resasc);

void gsl_integration_qk41 (const gsl_function * f, double a, double b,
			   double *result, double *abserr,
			   double *resabs, double *resasc);

void gsl_integration_qk51 (const gsl_function * f, double a, double b,
			   double *result, double *abserr,
			   double *resabs, double *resasc);

void gsl_integration_qk61 (const gsl_function * f, double a, double b,
			   double *result, double *abserr,
			   double *resabs, double *resasc);


/* The low-level integration rules in QUADPACK are identified by small
   integers (1-6). We'll use symbolic constants to refer to them.  */

enum
  {
    GSL_INTEG_GAUSS15 = 1,	/* 15 point Gauss-Kronrod rule */
    GSL_INTEG_GAUSS21 = 2,	/* 21 point Gauss-Kronrod rule */
    GSL_INTEG_GAUSS31 = 3,	/* 31 point Gauss-Kronrod rule */
    GSL_INTEG_GAUSS41 = 4,	/* 41 point Gauss-Kronrod rule */
    GSL_INTEG_GAUSS51 = 5,	/* 51 point Gauss-Kronrod rule */
    GSL_INTEG_GAUSS61 = 6	/* 61 point Gauss-Kronrod rule */
  };

int gsl_integration_qng (const gsl_function * f,
			 double a, double b,
			 double epsabs, double epsrel,
			 double *result, double *abserr,
			 size_t * neval);

int gsl_integration_qag (const gsl_function * f,
			 double a, double b,
			 double epsabs, double epsrel, size_t limit,
			 int key,
			 gsl_integration_workspace * workspace,
			 double *result, double *abserr);

int gsl_integration_qagi (gsl_function * f,
			  double epsabs, double epsrel, size_t limit,
			  gsl_integration_workspace * workspace,
			  double *result, double *abserr);

int gsl_integration_qagiu (gsl_function * f,
			   double a,
			   double epsabs, double epsrel, size_t limit,
			   gsl_integration_workspace * workspace,
			   double *result, double *abserr);

int gsl_integration_qagil (gsl_function * f,
			   double b,
			   double epsabs, double epsrel, size_t limit,
			   gsl_integration_workspace * workspace,
			   double *result, double *abserr);


int gsl_integration_qags (const gsl_function * f,
			  double a, double b,
			  double epsabs, double epsrel, size_t limit,
			  gsl_integration_workspace * workspace,
			  double *result, double *abserr);

int gsl_integration_qagp (const gsl_function * f,
			  double *pts, size_t npts,
			  double epsabs, double epsrel, size_t limit,
			  gsl_integration_workspace * workspace,
			  double *result, double *abserr);

int gsl_integration_qawc (gsl_function *f,
			  double a, double b, double c,
			  double epsabs, double epsrel, size_t limit,
			  gsl_integration_workspace * workspace,
			  double * result, double * abserr);

int gsl_integration_qaws (gsl_function * f,
			  const double a, const double b,
			  gsl_integration_qaws_table * t,
			  const double epsabs, const double epsrel,
			  const size_t limit,
			  gsl_integration_workspace * workspace,
			  double *result, double *abserr);

int gsl_integration_qawo (gsl_function * f,
			  const double a,
			  const double epsabs, const double epsrel,
			  const size_t limit,
			  gsl_integration_workspace * workspace,
			  gsl_integration_qawo_table * wf,
			  double *result, double *abserr);

int gsl_integration_qawf (gsl_function * f,
			  const double a,
			  const double epsabs,
			  const size_t limit,
			  gsl_integration_workspace * workspace,
			  gsl_integration_workspace * cycle_workspace,
			  gsl_integration_qawo_table * wf,
			  double *result, double *abserr);

#endif /* GSL_INTEGRATION_H */
