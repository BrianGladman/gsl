/* $Id$ */
#include <stdio.h>

#include "gsl_ran.h"
#include "gsl_randist.h"

main(int argc, char **argv)
{
    int i=1,N=10, *x;
    if (argc <= 1) {
	fprintf(stderr,"Usage: %s <N> [seed]\n",argv[0]);
	exit(0);
    } 
    if (argc > 1) N = atoi(argv[1]);
    if (argc > 2) {
	i = atoi(argv[2]);
	gsl_ran_seed(i);
    }
    fprintf(stderr,"N=%d\n",N);
    x = gsl_ran_shuffle(N,NULL);
    fprintf(stderr,"main: N=%d, x=%p\n",N,x);

    fprintf(stderr,"First shuffle:\n");
    for (i=0; i<N; ++i) {
	fprintf(stdout,"%d\n",x[i]);
    }
    fflush(stdout);
    /* shuffle again */
    fprintf(stderr,"Second shuffle:\n");
    gsl_ran_shuffle(N,x);
    for (i=0; i<N; ++i) {
	fprintf(stderr,"%d\n",x[i]);
    }
    fprintf(stderr,"Choose %d of them\n",N/2);
    gsl_ran_choose(N/2,N,x);
    for (i=0; i<N/2; ++i)
	fprintf(stderr,"%d\n",x[i]);
    fprintf(stderr,"And randomize those\n");
    gsl_ran_shuffle(N/2,x);
    for (i=0; i<N/2; ++i)
	fprintf(stderr,"%d\n",x[i]);
    /* free the memory */
    gsl_ran_shuffle(0,x);
    return 0;
}
