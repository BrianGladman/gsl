/* poly/init.c
 * 
 * Copyright (C) 2002 Gert Van den Eynde
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include <config.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_poly.h>

/*-*-*-*-*-*-*-*-*-*-*-* Allocators *-*-*-*-*-*-*-*-*-*-*-*/

gsl_poly * gsl_poly_alloc(const size_t size)
{
    gsl_poly * p = (gsl_poly *) malloc(sizeof(gsl_poly));

    if(p == 0) {
	GSL_ERROR_VAL("failed to allocate gsl_poly struct", GSL_ENOMEM, 0);
    }

    p->size = size;

    p->c = (double *) malloc(size * sizeof(double));

    if(p->c == 0) {
	free(p);
	GSL_ERROR_VAL("failed to allocate poly coefficients", GSL_ENOMEM, 0);
    }

    return p;
}

gsl_poly * gsl_poly_calloc(const size_t size)
{
    size_t k=0;
    gsl_poly * p = (gsl_poly *) malloc(sizeof(gsl_poly));

    if(p == 0) {
	GSL_ERROR_VAL("failed to allocate gsl_poly struct", GSL_ENOMEM, 0);
    }

    p->size = size;

    p->c = (double *) malloc(size * sizeof(double));

    if(p->c == 0) {
	free(p);
	GSL_ERROR_VAL("failed to allocate poly coefficients", GSL_ENOMEM, 0);
    }

    p->degree = 0;
    for (k=0;k<size;k++) p->c[k] = 0.0;

    return p;
}

int gsl_poly_set_Legendre(gsl_poly * p, unsigned int d, gsl_poly * w1, gsl_poly * w2, gsl_poly * w3) 
{
    gsl_poly * pk   = w1;
    gsl_poly * pkm1 = w2;
    gsl_poly * pkm2 = w3;
    gsl_poly * t;
    unsigned int k = 0, i = 0;

    if (p->size <= d) {
	GSL_ERROR_VAL("Size is not sufficient",GSL_EBADLEN,0);
    }

    gsl_poly_set_zero(p);

    if (d < 7) {
	switch (d) {
	    case 0:
		p->c[0]   = 1.0;
		p->degree = 0;
		break;
	    case 1:
		p->c[0]   = 0.0;
		p->c[1]   = 1.0;
		p->degree = 1;
		break;
	    case 2:
		p->c[0]   = -1.0/2.0;
		p->c[1]   = 0.0;
		p->c[2]   = 3.0/2.0;
		p->degree = 2;
		break;
	    case 3:
		p->c[0]   = 0.0;
		p->c[1]   = -3.0/2.0;
		p->c[2]   = 0.0;
		p->c[3]   = 5.0/2.0;
		p->degree = 3;
		break;
	    case 4:
		p->c[0]   = 3.0/8.0;
		p->c[1]   = 0.0;
		p->c[2]   = -30.0/8.0;
		p->c[3]   = 0.0;
		p->c[4]   = 35.0/8.0;
		p->degree = 4;
		break;
	    case 5:
		p->c[0]   = 0.0;
		p->c[1]   = 15.0/8.0;
		p->c[2]   = 0.0;
		p->c[3]   = -70.0/8.0;
		p->c[4]   = 0.0;
		p->c[5]   = 63.0/8.0;
		p->degree = 5;
		break;
	    case 6:
		p->c[0]   = -5.0/16.0;
		p->c[1]   = 0.0;
		p->c[2]   = 105.0/16.0;
		p->c[3]   = 0.0;
		p->c[4]   = -315.0/16.0;
		p->c[5]   = 0.0;
		p->c[6]   = 231.0/16.0;
		p->degree = 6;
		break;
	    default:
		GSL_ERROR_VAL("This should not happen",GSL_FAILURE,0);
		break;
	} 
    }
    else {

	if (pk->size < d+1)
	  {
	    GSL_ERROR_VAL("Size of first workspace poly not sufficient",GSL_EBADLEN,0);
	  }

	if (pkm1->size < d+1)
	  {
	    GSL_ERROR_VAL("Size of second workspace poly not sufficient",GSL_EBADLEN,0);
	  }

	if (pkm2->size < d+1)
	  {
	    GSL_ERROR_VAL("Size of third workspace poly not sufficient",GSL_EBADLEN,0);
	  }

	pkm2->c[0]   = 0.0;
	pkm2->c[1]   = 15.0/8.0;
	pkm2->c[2]   = 0.0;
	pkm2->c[3]   = -70.0/8.0;
	pkm2->c[4]   = 0.0;
	pkm2->c[5]   = 63.0/8.0;
	pkm2->degree = 5;

	pkm1->c[0]   = -5.0/16.0;
	pkm1->c[1]   = 0.0;
	pkm1->c[2]   = 105.0/16.0;
	pkm1->c[3]   = 0.0;
	pkm1->c[4]   = -315.0/16.0;
	pkm1->c[5]   = 0.0;
	pkm1->c[6]   = 231.0/16.0;
	pkm1->degree = 6;

	for (k = 7; k <= d; k++) 
	  {
	    gsl_poly_set_zero(pk);

	    for (i = 1; i <= k; i++) 
	      {
		pk->c[i] = pkm1->c[i-1];
	      }
	    gsl_poly_set_degree(pk, 100.0*GSL_DBL_EPSILON);

	    gsl_poly_scale(pk, (2.0 * k - 1.0) / k);

	    gsl_poly_scale(pkm2, (k - 1.0) / k);

	    gsl_poly_sub(pk,pkm2);

	    t    = pkm2;
	    pkm2 = pkm1;
	    pkm1 = pk;
	    pk   = t;

	  }

	gsl_poly_memcpy(p,pkm1);

    }
    return GSL_SUCCESS;
}

void gsl_poly_free(gsl_poly * p)
{
    free(p->c);
    free(p);
}


gsl_poly_sturm * gsl_poly_sturm_alloc(const size_t size)
{
    int i = 0, k = 0;
    gsl_poly_sturm * ss = (gsl_poly_sturm *) malloc(sizeof(gsl_poly_sturm));

    if (ss == 0) {
	GSL_ERROR_VAL("failed to allocate gsl_poly_sturm struct", GSL_ENOMEM, 0);
    }

    ss->sturmseq = malloc(size*sizeof(gsl_poly *));

    if (ss->sturmseq == 0) {
	free(ss);
	GSL_ERROR_VAL("failed to allocate gsl_poly_sturm data array", GSL_ENOMEM, 0);
    }

    for (i = 0; i < size; i++) {
	ss->sturmseq[i] = gsl_poly_alloc(size);

	if (ss->sturmseq[i] == 0) {
	    for (k = 0; k < i; k++) {
		gsl_poly_free(ss->sturmseq[k]);
	    }
	    free(ss);
	    GSL_ERROR_VAL("failed to allocate gsl_poly_sturm data array", GSL_ENOMEM, 0);

	}
    }
    ss->size = size;

    return ss;
} 

gsl_poly_sturm * gsl_poly_sturm_calloc(const size_t size)
{
    int i = 0, k = 0;
    gsl_poly_sturm * ss = (gsl_poly_sturm *) malloc(sizeof(gsl_poly_sturm));

    if (ss == 0) {
	GSL_ERROR_VAL("failed to allocate gsl_poly_sturm struct", GSL_ENOMEM, 0);
    }

    ss->sturmseq = malloc(size*sizeof(gsl_poly *));

    if (ss->sturmseq == 0) {
	free(ss);
	GSL_ERROR_VAL("failed to allocate gsl_poly_sturm data array", GSL_ENOMEM, 0);
    }

    for (i = 0; i < size; i++) {
	ss->sturmseq[i] = gsl_poly_calloc(size);

	if (ss->sturmseq[i] == 0) {
	    for (k = 0; k < i; k++) {
		gsl_poly_free(ss->sturmseq[k]);
	    }
	    free(ss);
	    GSL_ERROR_VAL("failed to allocate gsl_poly_sturm data array", GSL_ENOMEM, 0);
	}
    }
    ss->size = size;

    return ss;
} 

void gsl_poly_sturm_free(gsl_poly_sturm * ss) 
{
    int i = 0;
    int size = ss->size;

    for (i = 0; i < size; i++) {
	gsl_poly_free(ss->sturmseq[i]);
    }

    free(ss->sturmseq);

    free(ss);

} 


