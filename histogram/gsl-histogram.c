#include <stdio.h>
#include <stdlib.h>
#include <gsl_histogram.h>

int
main (int argc, char **argv)
{
  double a = 0.0, b = 1.0 ;
  size_t n = 10;

  if (argc != 4)
    {
      printf ("Usage: gsl-histogram xmin xmax n\n");
      printf ("Computes a histogram of the data on stdin using n bins from xmin to xmax\n");
      exit (0);
    }

  a = atof (argv[1]);
  b = atof (argv[2]);
  n = atoi (argv[3]);

  {
    gsl_histogram * h = gsl_histogram_calloc_uniform (n, a, b) ;
    int status ;
    
    do {
      double x ;
      status = fscanf(stdin, "%lg", &x) ;
      
      gsl_histogram_increment (h, x) ;

    } while (status == 1) ;

    gsl_histogram_fprintf (stdout, h, "%g") ;
  }
  
  return 0 ;
}
