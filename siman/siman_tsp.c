#include <math.h>

#include <gsl_ran.h>
#include <gsl_siman.h>
#include <stdio.h>

/* set up parameters for this simulated annealing run */
#define N_TRIES 200		/* how many points do we try before stepping */
#define ITERS_FIXED_T 1000	/* how many iterations for each T? */
#define STEP_SIZE 1.0		/* max step size in random walk */
#define K 1.0			/* Boltzmann constant */
#define T_INITIAL 80.1		/* initial temperature */
#define MU_T 1.003		/* damping factor for temperature */
#define T_MIN 1.0e-2

gsl_siman_params_t params = {N_TRIES, ITERS_FIXED_T, STEP_SIZE,
			     K, T_INITIAL, MU_T, T_MIN};

#define N_CITIES 10

#define CITY_NAME_LEN 50
struct s_tsp_city {
  char name[CITY_NAME_LEN];
  double x, y;			/* coordinates */
};
typedef struct s_tsp_city Stsp_city;

Stsp_city cities[N_CITIES] = {{"Santa_Fe",    10.2,	35.4},
			      {"Albuquerque", 12.1,	31.4},
			      {"Clovis",       6.9,	30.4},
			      {"Dallas",       4.2,	29.4},
			      {"Alamagordo",  12.1,	26.1},
			      {"Los_Alamos",  12.9,	36.0},
			      {"Tesuque",     10.2,	35.7},
			      {"Las_Cruces",  12.1,	18.4},
			      {"Phoenix",     30.1,	22.4},
			      {"Tuscon",      31.2,	20.0}};

/* distance between two cities */
double city_distance(Stsp_city c1, Stsp_city c2)
{
  double dx = c2.x - c1.x;
  double dy = c2.y - c1.y;
  return sqrt(dx*dx + dy*dy);
}

/* energy for the travelling salesman problem */
double Etsp(void *xp)
{
  /* an array of N_CITIES integers describing the order */
  int *route = (int *) xp;
  double E = 0;
  int i;

  for (i = 0; i < N_CITIES; ++i) {
/*     E += city_distance_matrix[route[i]][route[(i + 1) % N_CITIES]]; */
    E += city_distance(cities[route[i]], cities[route[(i + 1) % N_CITIES]]);
  }

  return E;
}

double Mtsp(void *xp, void *yp)
{
  int *route1 = (int *) xp, *route2 = (int *) yp;
  double distance = 0;
  int i;

  for (i = 0; i < N_CITIES; ++i) {
    distance += ((route1[i] == route2[i]) ? 0 : 1);
  }

  return distance;
}

/* take a step through the TSP space */
void Stsp(void *xp, double step_size)
{
  int x1, x2, dummy;
  int *route = (int *) xp;

  /* pick the two cities to swap in the matrix */
  x1 = gsl_ran_random() % N_CITIES;
  do {
    x2 = gsl_ran_random() % N_CITIES;
  } while (x2 == x1);

  dummy = route[x1];
  route[x1] = route[x2];
  route[x2] = dummy;
}

void Ptsp(void *xp)
{
  int i;
  int *route = (int *) xp;
  printf("  [");
  for (i = 0; i < N_CITIES; ++i) {
    printf(" %d ", route[i]);
  }
  printf("]  ");
}

int main(int argc, char *argv[])
{
  int x_initial[N_CITIES];
  int i;

  /* set up a trivial initial route */
  printf("# initial order of cities:\n");
  for (i = 0; i < N_CITIES; ++i) {
    printf("# %s\n", cities[i].name);
    x_initial[i] = i;
  }

/*   gsl_siman_solve(x_initial, Etsp, Stsp, Mtsp, NULL, */
/* 		  N_CITIES*sizeof(int), params); */

  gsl_siman_solve(x_initial, Etsp, Stsp, Mtsp, Ptsp,
		  N_CITIES*sizeof(int), params);

  printf("# final order of cities:\n");
  for (i = 0; i < N_CITIES; ++i) {
    printf("# %s\n", cities[x_initial[i]].name);
  }

  return 0;
}
