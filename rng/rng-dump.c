#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl_rng.h>

int
main (int argc, char **argv)
{
  int i, j ;
  unsigned long buffer[1024] ;
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


  for (i = 0; i < 3000 ; i++)
    {
      int status ;

      for (j = 0; j < 1024; j++)
	{
	  buffer[j] = gsl_rng_get (r) ;
	}
      
      status = fwrite(buffer, sizeof(unsigned long), 1024, stdout) ;

      if (status != 1024) 
	{
	  perror("fwrite") ;
	  exit(EXIT_FAILURE) ;
	}
    }

  return 0;
}
