/* poly/sturm.c
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

int gsl_poly_sturm_build(gsl_poly_sturm * ss, const gsl_poly * p, gsl_poly * w)
{
    int d = p->degree;
    int index = 0;
    double alpha = 0.0;

    gsl_poly * q = w;

    if (q->size <= p->degree) {
	GSL_ERROR("Size of workspace poly not sufficient",GSL_EBADLEN);
    }


    alpha = gsl_poly_get(p,d); 

    if(ss->size < p->degree) {
	free(q);
	GSL_ERROR("Sturm sequence not large enough",GSL_EINVAL);
    }

    gsl_poly_memcpy(ss->sturmseq[0],p);

    gsl_poly_scale(ss->sturmseq[0],alpha);

    gsl_poly_diff(ss->sturmseq[0],ss->sturmseq[1]);

    alpha = fabs(gsl_poly_get(ss->sturmseq[1],d-1));
    gsl_poly_scale(ss->sturmseq[1],1.0/alpha);

    index = 2;

    while (ss->sturmseq[index-1]->degree > 0) {
	gsl_poly_div(ss->sturmseq[index-2],ss->sturmseq[index-1],q,ss->sturmseq[index]);
	d = ss->sturmseq[index]->degree;
	alpha = fabs(gsl_poly_get(ss->sturmseq[index],d));
	if (ss->sturmseq[index]->degree > 0) {
	    gsl_poly_scale(ss->sturmseq[index],-1.0/alpha);
	}
	else {
	    gsl_poly_scale(ss->sturmseq[index],-1.0);
	} 
	index++;
    }

    ss->index = index - 1;
    return GSL_SUCCESS;
}

int gsl_poly_sturm_changes(const gsl_poly_sturm * ss, double a) 
{
    int changes = 0;
    double u = 0.0, v = 0.0;
    int i = 0;
    int index = ss->index;

    v = gsl_poly_eval2(ss->sturmseq[0],a);

    for (i = 1; i <= index; i++) {
	u = gsl_poly_eval2(ss->sturmseq[i],a);
	if (v == 0.0 || v * u < 0.0) {
	    changes++;
	}
	v = u;
    }


    return changes;

}

int gsl_poly_sturm_numroots(const gsl_poly_sturm * ss,double a, double b)
{
    int nr = gsl_poly_sturm_changes(ss,a) - gsl_poly_sturm_changes(ss,b);

    return nr;
}

void gsl_poly_sturm_dump(const gsl_poly_sturm * ss)
{

    int i = 0,j = 0;
    gsl_poly * p;

    for (i = 0;i<=ss->index;i++) {
	p = ss->sturmseq[i];
	for (j=p->degree;j>=0;j--) {
	    printf("%f ",p->c[j]);
	}
	printf("\n");
    }
}
