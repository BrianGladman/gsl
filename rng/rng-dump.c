#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl_rng.h>

int
main (int argc, char **argv)
{
  int i, j ;
  char buffer[1024 * 4] ;
  gsl_rng * r ;

  gsl_rng_env_setup () ;

  r = gsl_rng_alloc (gsl_rng_default) ;

  if (argc != 1)
    {
      printf ("Usage: diehard\n");
      printf ("Output 3 million numbers in binary format," 
	      "suitable for testing with DIEHARD\n");
      exit (0);
    }

  argv = 0 ; /* prevent warning about unused argument */

  for (i = 0; i < 3000 ; i++)
    {
      int status ;

      for (j = 0; j < 1024; j++)
	{
	  unsigned long int u = gsl_rng_get (r) ;
	  buffer[4 * j + 0] = u & 0xFF ;
	  u >>= 8;
	  buffer[4 * j + 1] = u & 0xFF ;
	  u >>= 8;
	  buffer[4 * j + 2] = u & 0xFF ;
	  u >>= 8;
	  buffer[4 * j + 3] = u & 0xFF ;

	}
      
      status = fwrite(buffer, 4 * sizeof(char), 1024, stdout) ;

      if (status != 1024) 
	{
	  perror("fwrite") ;
	  exit(EXIT_FAILURE) ;
	}
    }

  return 0;
}
