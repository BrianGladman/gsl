#include "qpsrt.c"

static inline
void initialise (gsl_integration_workspace * workspace, double a, double b);

static inline
void set_initial_result (gsl_integration_workspace * workspace, 
			 double result, double error);

static inline
void update (gsl_integration_workspace * workspace,
		 double a1, double b1, double area1, double error1,
		 double a2, double b2, double area2, double error2);

static inline void
retrieve (const gsl_integration_workspace * workspace, 
	  double * a, double * b, double * r, double * e);

static double
sum_results (const gsl_integration_workspace * workspace);


static inline
void initialise (gsl_integration_workspace * workspace, double a, double b)
{
  workspace->size = 0;
  workspace->nrmax = 0;
  workspace->i = 0;
  workspace->alist[0] = a;
  workspace->blist[0] = b;
  workspace->rlist[0] = 0.0;
  workspace->elist[0] = 0.0;
  workspace->order[0] = 0;
  workspace->level[0] = 0;

  workspace->maximum_level = 0;
}

static inline
void set_initial_result (gsl_integration_workspace * workspace, 
			 double result, double error)
{
  workspace->size = 1;
  workspace->rlist[0] = result;
  workspace->elist[0] = error;
}

static inline
void append_interval (gsl_integration_workspace * workspace,
		      double a1, double b1, double area1, double error1)
{
  const size_t i_new = workspace->size ;

  workspace->alist[i_new] = a1;
  workspace->blist[i_new] = b1;
  workspace->rlist[i_new] = area1;
  workspace->elist[i_new] = error1;
  workspace->order[i_new] = i_new;
  workspace->level[i_new] = 0;

  workspace->size++;
}


static inline
void update (gsl_integration_workspace * workspace,
	     double a1, double b1, double area1, double error1,
	     double a2, double b2, double area2, double error2)
{
  double * alist = workspace->alist ;
  double * blist = workspace->blist ;
  double * rlist = workspace->rlist ;
  double * elist = workspace->elist ;
  size_t * level = workspace->level ;

  const size_t i_max = workspace->i ;
  const size_t i_new = workspace->size ;

  const size_t new_level = workspace->level[i_max] + 1;

  /* append the newly-created intervals to the list */
  
  if (error2 > error1)
    {
      alist[i_max] = a2;	/* blist[maxerr] is already == b2 */
      rlist[i_max] = area2;
      elist[i_max] = error2;
      level[i_max] = new_level;
      
      alist[i_new] = a1;
      blist[i_new] = b1;
      rlist[i_new] = area1;
      elist[i_new] = error1;
      level[i_new] = new_level;
    }
  else
    {
      blist[i_max] = b1;	/* alist[maxerr] is already == a1 */
      rlist[i_max] = area1;
      elist[i_max] = error1;
      level[i_max] = new_level;
      
      alist[i_new] = a2;
      blist[i_new] = b2;
      rlist[i_new] = area2;
      elist[i_new] = error2;
      level[i_new] = new_level;
    }
  
  workspace->size++;

  if (new_level > workspace->maximum_level)
    {
      workspace->maximum_level = new_level;
    }

  qpsrt (workspace) ;
}

static double
sum_results (const gsl_integration_workspace * workspace)
{
  const double * const rlist = workspace->rlist ;
  const size_t n = workspace->size;

  size_t k;
  double result_sum = 0;

  for (k = 0; k < n; k++)
    {
      result_sum += rlist[k];
    }
  
  return result_sum;
}

static inline void
retrieve (const gsl_integration_workspace * workspace, 
	  double * a, double * b, double * r, double * e)
{
  const size_t i = workspace->i;
  double * alist = workspace->alist;
  double * blist = workspace->blist;
  double * rlist = workspace->rlist;
  double * elist = workspace->elist;

  *a = alist[i] ;
  *b = blist[i] ;
  *r = rlist[i] ;
  *e = elist[i] ;
}

static inline void
reset_nrmax (gsl_integration_workspace * workspace);


static inline void
reset_nrmax (gsl_integration_workspace * workspace)
{
  workspace->nrmax = 0;
  workspace->i = workspace->order[0] ;
}

static void
sort_results (gsl_integration_workspace * workspace);

static void
sort_results (gsl_integration_workspace * workspace)
{
  size_t i;
  
  double * elist = workspace->elist ;
  size_t * order = workspace->order ;

  size_t nint = workspace->size;

  for (i = 0; i < nint; i++)
    {
      size_t i1 = order[i];
      double e1 = elist[i1];
      size_t i_max = i1;
      size_t j;

      for (j = i + 1; j < nint; j++)
	{
	  size_t i2 = order[j];
	  double e2 = elist[i2];

	  if (e2 >= e1)
	    {
	      i_max = i2;
	      e1 = e2;
	    }
	}

      if (i_max != i1)
	{
	  order[i] = order[i_max];
	  order[i_max] = i1;
	}
    }

  workspace->i = order[0] ;
}

static inline int
subinterval_too_small (double a1, double a2, double b2);


static inline int
subinterval_too_small (double a1, double a2, double b2)
{
  const double e = GSL_DBL_EPSILON;
  const double u = GSL_DBL_MIN;

  double tmp = (1 + 100 * e) * (fabs (a2) + 1000 * u);

  int status = fabs (a1) <= tmp && fabs (b2) <= tmp;

  return status;
}

/* Compare the integral of f(x) with the integral of |f(x)|
   to determine if f(x) covers both positive and negative values */

static inline int
test_positivity (double result, double resabs);

static inline int
test_positivity (double result, double resabs)
{
  int status = (fabs (result) >= (1 - 50 * GSL_DBL_EPSILON) * resabs);

  return status;
}
