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

#define CITY_NAME_LEN 50
struct s_tsp_city {
  char name[CITY_NAME_LEN];
  double lat, longitude;	/* coordinates */
};
typedef struct s_tsp_city Stsp_city;

/* in this table, latitude and longitude are obtained from the US
   Census Bureau, at http://www.census.gov/cgi-bin/gazetteer */
Stsp_city cities[] = {{"Santa_Fe",    35.68,   105.95},
		      {"Albuquerque", 35.12,   106.62},
		      {"Clovis",      34.41,   103.20},
		      {"Dallas",      32.79,    96.77},
		      {"Grants",      35.15,   107.84},
		      {"Los Alamos",  35.89,   106.28},
		      {"Tesuque",     35.77,   105.92},
		      {"Las Cruces",  32.34,   106.76},
		      {"Phoenix",     33.54,   112.07},
		      {"Durango",     37.29,   107.87},
		      {"Cortez",      37.35,   108.58},
		      {"Gallup",      35.52,   108.74}};
/*  			      {"Tuscon",      }}; */

#define N_CITIES (sizeof(cities)/sizeof(Stsp_city))


/* distance between two cities */
double city_distance(Stsp_city c1, Stsp_city c2)
{
  const earth_radius = 7;	/* 7KM approximately */
  double x1 = earth_radius*cos(2*c1.lat)*cos(c1.longitude);
  double x2 = earth_radius*cos(2*c2.lat)*cos(c2.longitude);

  double y1 = earth_radius*cos(2*c1.lat)*sin(c1.longitude);
  double y2 = earth_radius*cos(2*c2.lat)*sin(c2.longitude);

  double z1 = earth_radius*sin(2*c1.lat);
  double z2 = earth_radius*sin(2*c2.lat);

  double dx = x2 - x1;
  double dy = y2 - y1;
  double dz = z2 - z1;

  return sqrt(dx*dx + dy*dy + dz*dz);
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
