

/* $Id$ */
/* Driver routine for the uniform random number generators */
#include <stdio.h>
#include <stdlib.h>
#include "gsl_ran.h"
#include "gsl_randist.h"
int
main (int argc, char **argv)
{
  int i, n = 1;
  int randseed = 17;
  if (argc == 1)
    {
      printf ("Usage: %s <n> [seed]\n", argv[0]);
      printf ("          Writes <n> random numbers\n");
      printf ("          normally distributed mean=0, var=1,\n");
      printf ("          optionally using <seed>\n");
      exit (0);
    }
  if (argc > 1)
    n = atoi (argv[1]);
  if (argc > 2)
    {
      randseed = atoi (argv[2]);
      gsl_ran_seed (randseed);
    }

  for (i = 0; i < n; ++i)
    {
      printf ("%.10f\n", gsl_ran_gaussian ());
    }
  return 0;
}
