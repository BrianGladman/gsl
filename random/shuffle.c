/* $Id$ */
/* Randomly permute (shuffle) N indices */
/* Supply an integer array x[N], and on return, it will
 * be filled with indices 0...N-1 in random order.
 * The algorithm is from Knuth, SemiNumerical Algorithms, v2, p139 
 */
#include <stdlib.h>
#include <stdio.h>		/* defines NULL */
#include <math.h>		/* defines floor() */
#include "gsl_ran.h"

int *
gsl_ran_shuffle(int N, int *x)
{
    int i,k,tmp;

    /* First, do a bunch of memory allocation stuff */
    if (N<0) return NULL;
    if (x==NULL && N>0) {
	x = (int *)calloc(N,sizeof(int));
	if (x==NULL) return NULL;
	for (i=0; i<N; ++i)	
	    x[i]=i;
    }
    if (x != NULL && N==0) {
	cfree((char *)x);
	return NULL;
    }
    /* Now here's the algorithm, more or less transcribed
     * from Knuth, who cites Moses and Oakford, and Durstenfeld */

    for (i=N-1; i>=0; --i) {
	k = floor(i*gsl_ran_uniform());
	tmp = x[k]; x[k]=x[i]; x[i]=tmp;
    }
    return x;
}
int *
gsl_ran_choose(int K, int N, int *x)
{
    int n,k;
    /* Choose K out of N items */
    /* return an array x[] of the indices of the N items */
    /* these items will be in sorted order -- you can use
     * shuffle() to randomize them if you wish */

    /* First, do a bunch of memory allocation stuff */
    if (N<K || K<0) {
	if (x != NULL)
	    cfree((char *)x);
	return NULL;
    }
    if (x==NULL && K>0) {
	x = (int *)calloc(K,sizeof(int));
	if (x==NULL) 
	    return NULL;
    }
    /* Here is the guts of the algorithm: three lines!! */
    for (n=0, k=0; n<N && k<K; ++n) {
	if ((N-n)*gsl_ran_uniform() < K-k) {
	    x[k++] = n;
	}
    }
    return x;
}
	
#ifdef MAIN
main(int argc, char **argv)
{
    int i=1,N=10, *x;
    if (argc <= 1) {
	fprintf(stderr,"Usage: %s <N>\n",argv[0]);
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
#endif
    

