#ifndef __SIMAN_H_
#define __SIMAN_H_

/* this structure contains all the information needed to structure
 * the search, beyond the energy function, the step function and
 * the initial guess.
 */
struct s_siman_params {
  int n_tries;		/* how many points to try for each step */
  int iters_fixed_T;	/* how many iterations at each temperature? */
  double step_size;	/* max step size in the random walk */
  /* the following parameters are for the Boltzmann distribution */
  double k, t_initial, mu_t, t_min;
};

typedef struct s_siman_params Ssiman_params;



struct CA_rule {
  char *str;
  int k, n, length;
};

union u_Element {
  double D1;
  double D2[2];
  double D3[3];
  struct CA_rule rule;
};

typedef union u_Element Element;

/* prototype for the workhorse function */
void gsl_siman_solve(Element *x0_p, double (*Efunc)(Element x),
		     void (*take_step)(Element *x_p, double step_size),
		     double distance(Element x, Element y),
		     void print_position(Element x),
		     Ssiman_params params);


/* 1-dimensional energy and stepping functions */
double test_E_1D(Element x);
void test_step_1D(Element *x_p, double step_size);
double distance_1D(Element x, Element y);
void print_pos_1D(Element x);

/* 2-dimensional energy and stepping functions */
double test_E_2D(Element x);
void test_step_2D(Element *x_p, double step_size);
double distance_2D(Element x, Element y);
void print_pos_2D(Element x);

/* 3-dimensional energy and stepping functions */
double test_E_3D(Element x);
void test_step_3D(Element *x_p, double step_size);
double distance_3D(Element x, Element y);
void print_pos_3D(Element x);

void siman_solve(Element *x0_p, double (*Efunc)(Element x),
		 void (*take_step)(Element *x_p, double step_size),
		 double distance(Element x, Element y),
		 void print_position(Element x),
		 Ssiman_params params);

#endif /* __SIMAN_H */
