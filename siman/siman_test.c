#include <math.h>

#include <gsl_ran.h>
#include <gsl_siman.h>
#include <stdio.h>

/* set up parameters for this simulated annealing run */
#define N_TRIES 200		/* how many points do we try before stepping */
#define ITERS_FIXED_T 1000	/* how many iterations for each T? */
#define STEP_SIZE 1.0		/* max step size in random walk */
#define K 1.0			/* Boltzmann constant */
#define T_INITIAL 0.008		/* initial temperature */
#define MU_T 1.003		/* damping factor for temperature */
#define T_MIN 2.0e-6

gsl_siman_params_t params = {N_TRIES, ITERS_FIXED_T, STEP_SIZE,
			     K, T_INITIAL, MU_T, T_MIN};

double square(double x);

double test_E_1D(Element x)
{
  double val = x.D1;
/*return sin(sin(val*val) - cos(val)) + cos(sin(val) + sin(val)*sin(val));*/
  return exp(-square(val-1))*sin(8*val);
/*  return 1.0/(square(x-1.0) + 1.0)*sin(8.0*x); */
/*  return 1.0/(square(x-1.0) + 1.0)*sin(8.0*x) + x*x/100.0; */
}

double square(double x)
{
  return x*x;
}

/* takes a step for the test function; max distance: step_size.
 * the new point is put in x_p and returned.
 */
void test_step_1D(Element *x_p, double step_size)
{
  double r;
  double old_x = x_p->D1;
  double new_x;

  r = gsl_ran_uniform();
  new_x = r;
  new_x = new_x*2*params.step_size;
  new_x = new_x - params.step_size + old_x;

  x_p->D1 = new_x;
}


/* simple routine to print out a position value */
void print_pos_1D(Element x)
{
  printf("%12g", x.D1);
}

/* a metric for the 2D space */
double distance_1D(Element x, Element y)
{
  return fabs(y.D1 - x.D1);
}


/* a 2-D function to be minimized */
double test_E_2D(Element x)
{
  double old_x = x.D2[0], old_y = x.D2[1];

  return exp(-square(old_x-1) - square(old_y - 0.8))*sin(8*old_x + 8 * old_y);
}

/* takes a step for the test function; max distance: step_size.  the
   new point is put in x_p and returned. */
void test_step_2D(Element *x_p, double step_size)
{
  double r;
  double old_x =  x_p->D2[0], old_y = x_p->D2[1], new_x, new_y;

  r = gsl_ran_uniform();
  new_x = r;
  new_x = new_x*2*params.step_size;
  new_x = new_x - params.step_size + old_x;

  r = gsl_ran_uniform();
  new_y = r;
  new_y = new_y*2*params.step_size;
  new_y = new_y - params.step_size + old_y;

  x_p->D2[0] = new_x;
  x_p->D2[1] = new_y;
}

/* simple routine to print out a position value */
void print_pos_2D(Element x)
{
  printf("%g::%g", x.D2[0], x.D2[1]);
}

/* a metric for the 2D space */
double distance_2D(Element x, Element y)
{
  return sqrt(square(y.D2[0]-x.D2[0]) + square(y.D2[1]-x.D2[1]));
}


/**********************************************/
/************ 3-dimensional search ************/
/**********************************************/

/* a 3-D function to be minimized */
double test_E_3D(Element x)
{
  return exp(-square(x.D3[0]-1) - square(x.D3[1] - 0.8)
    - square(x.D3[2] - 0.8)) * sin(8*x.D3[0] + 8*x.D3[1] + 8*x.D3[2])
      + (square(x.D3[0]) + square(x.D3[1]) + square(x.D3[2]))/10000.0;
}

/* takes a step for the test function; max distance: step_size.
 * the new point is put in x_p and returned.
 */
void test_step_3D(Element *x_p, double step_size)
{
  double r;
  double old_x = x_p->D3[0], old_y = x_p->D3[1], old_z = x_p->D3[2];
  double new_x, new_y, new_z;

  r = gsl_ran_uniform();
  new_x = r;
  new_x = new_x*2*params.step_size;
  new_x = new_x - params.step_size + old_x;

  r = gsl_ran_uniform();
  new_y = r;
  new_y = new_y*2*params.step_size;
  new_y = new_y - params.step_size + old_y;

  r = gsl_ran_uniform();
  new_z = r;
  new_z = new_z*2*params.step_size;
  new_z = new_z - params.step_size + old_z;

  x_p->D3[0] = new_x;
  x_p->D3[1] = new_y;
  x_p->D3[2] = new_z;
}

/* simple routine to print out a position value */
void print_pos_3D(Element x)
{
  printf("%g::%g::%g", x.D3[0], x.D3[1], x.D3[2]);
}

/* a metric for the 2D space */
double distance_3D(Element x, Element y)
{
  return sqrt(square(y.D3[0]-x.D3[0]) + square(y.D3[1]-x.D3[1])
	      + square(y.D3[2]-x.D3[2]));
}

/* now some functions to test in one dimension */
double E1(void *xp)
{
  double x = * ((double *) xp);

/*   return exp(-square(x-1))*sin(8*x); */
  return exp(-square(x-1))*sin(8*x) - exp(-square(x-1000))*0.89;
}

double M1(void *xp, void *yp)
{
  double x = *((double *) xp);
  double y = *((double *) yp);

  return fabs(x - y);
}

void S1(void *xp, double step_size)
{
  double r;
  double old_x = *((double *) xp);
  double new_x;

  r = gsl_ran_uniform();
  new_x = r*2*step_size - step_size + old_x;
/*   new_x = new_x*2*step_size; */
/*   new_x = new_x - step_size + old_x; */

  memcpy(xp, &new_x, sizeof(new_x));
}

void P1(void *xp)
{
  printf(" %12g ", *((double *) xp));
}

int main(int argc, char *argv[])
{
  Element x0;			/* initial guess for search */

/*   double x_initial = 2.5; */
  double x_initial = -10.0;

  gsl_siman_solve(&x_initial, E1, S1, M1, P1, sizeof(double), params);

  return 0;

  if (argc != 2) {
    fprintf(stderr, "usage: %s [D1, D2, D3, CA]\n", argv[0]);
    return 1;
  }

  printf("#testing the simulated annealing routines\n");
  if (strcmp(argv[1], "D1") == 0) {
    x0.D1 = 12.0;
    printf("#one dimensional problem, x0 = %f\n", x0.D1);
/*     gsl_siman_Usolve(&x0, test_E_1D, test_step_1D, distance_1D, */
/* 		    print_pos_1D, params); */
    return 0;
  }

  if (strcmp(argv[1], "D2") == 0) {
    x0.D2[0] = 12.0;
    x0.D2[1] = 5.5;
    printf("#two dimensional problem, (x0,y0) = (%f,%f)\n",
	   x0.D2[0], x0.D2[1]);
/*     gsl_siman_Usolve(&x0, test_E_2D, test_step_2D, distance_2D, */
/* 		    print_pos_2D, params); */
    return 0;
  }

  if (strcmp(argv[1], "D3") == 0) {
    x0.D3[0] = 12.2;
    x0.D3[1] = 5.5;
    x0.D3[2] = -15.5;
    printf("#three dimensional problem, (x0,y0,z0) = (%f,%f,%f)\n",
	   x0.D3[0], x0.D3[1], x0.D3[2]);
/*     gsl_siman_Usolve(&x0, test_E_3D, test_step_3D, distance_3D, */
/* 		    print_pos_3D, params); */
  }
/*
  x0.D2[0] = 12.2;
  x0.D2[1] = 5.5;

  gsl_siman_solve(&x0, test_E_2D, test_step_2D, distance_2D, print_pos_2D, params);
*/
/*
  x0.D3[0] = 12.2;
  x0.D3[1] = 5.5;
  x0.D3[2] = -15.5;

  gsl_siman_solve(&x0, test_E_3D, test_step_3D, distance_3D, print_pos_3D, params);
  */

  return 0;
}
