/* poly/oper.c
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

extern int gsl_check_range ; /* defined in vector/vector.c */

void gsl_poly_dump (gsl_poly * p)
{
    int i;
    printf("%d\n",p->degree);
    for (i = p->degree; i >= 0; i--) 
      {
	printf("%e\n",p->c[i]);
      }
    printf("\n");
}

double gsl_poly_get(const gsl_poly * p, const size_t i)
{
    if (gsl_check_range)
      {
	if (i > p->degree)
	  {
	    GSL_ERROR_VAL("Index out of range", GSL_EINVAL, 0.0);
	  }
      }
    return p->c[i];
}

void gsl_poly_set(gsl_poly * p, const size_t i, double a)
{
    if (gsl_check_range)
      {
	if (i > p->size - 1)
	  {
	    GSL_ERROR_VOID("Index out of range", GSL_EINVAL);
	  }
      }
    p->c[i] = a;
    if (i > p->degree) p->degree = i;
}

int gsl_poly_set_zero(gsl_poly * p)
{ 
    size_t i;

    for (i = 0; i < p->size; i++)
      {
	p->c[i] = 0.0;
      }
    p->degree = 0;
    return GSL_SUCCESS;
}

int gsl_poly_set_all(gsl_poly * p, size_t d, double x)
{
    size_t i;
    if (p->size <= d) 
      { 
	GSL_ERROR("size of destination poly too small",GSL_EBADLEN);
      }	
    else
      {    
	p->degree = d;
	for (i = 0; i <= d; i++)
	  {
	    p->c[i] = x;
	  }
      }
    return GSL_SUCCESS;
}

int gsl_poly_consistent(const gsl_poly * p, double tol) 
{
    int i = 0, d = 0;
    int status = 0;

    for (i = p->size - 1; i >= 0; i--)
      {
	if (fabs(p->c[i]) >= tol)
	  {
	    d = i;
	    break;
	  }
      }
    if ( d == p->degree) 
      {
	status = GSL_SUCCESS;
      }
    else
      {
	status = GSL_FAILURE;
      }
    return status;
}

int gsl_poly_scale(gsl_poly * p, double a)
{
    size_t d = p->degree;
    size_t i;

    for (i = 0; i <= d; i++)
      {
	p->c[i] *= a;
      }

    return GSL_SUCCESS;
}

int gsl_poly_memcpy(gsl_poly * dest, const gsl_poly * src)
{
    size_t src_size = src->size;
    size_t dest_size = dest->size;
    size_t j;


    if (src_size > dest_size)
      {
	GSL_ERROR("size of destination poly too small",GSL_EBADLEN);
      }	

    for (j = 0; j < src_size; j++)
      {
	dest->c[j] = src->c[j];
      }
    dest->degree = src->degree;

    return GSL_SUCCESS;
}

double gsl_poly_eval2(gsl_poly * p, double x)
{
    return gsl_poly_eval(p->c,p->degree+1,x);
}

int gsl_poly_get_degree(const gsl_poly * p)
{
    return p->degree;
}

int gsl_poly_set_degree(gsl_poly * p, double tol) 
{
    int i = 0;
    size_t d = 0;

    for (i = p->size - 1; i >= 0; i--)
      {
	if (fabs(p->c[i]) < tol)
	  {
	    p->c[i] = 0.0;
	  }
	else
	  {
	    d = i;
	    break;
	  }
      }
    p->degree = d;

    return (GSL_SUCCESS);
}


int gsl_poly_add(gsl_poly * p1, const gsl_poly * p2)
{
    size_t i = 0;

    if (p1->size <= p2->degree)
      {
	GSL_ERROR("size of destination poly not sufficient", GSL_EBADLEN);
      }
    else
      {
	for (i = 0; i <= p2->degree; i++)
	  {
	    p1->c[i] += p2->c[i];
	  }
	p1->degree = GSL_MAX(p1->degree, p2->degree);
	return GSL_SUCCESS;
      }
}

int gsl_poly_sub(gsl_poly * p1, const gsl_poly * p2)
{
    size_t i = 0;

    if (p1->size <= p2->degree)
      {
	GSL_ERROR("size of destination poly not sufficient", GSL_EBADLEN);
      }
    else
      {
	for (i = 0; i <= p2->degree; i++)
	  {
	    p1->c[i] -= p2->c[i];
	  }
	p1->degree = GSL_MAX(p1->degree, p2->degree);
	return GSL_SUCCESS;
      }
}

int gsl_poly_mul(gsl_poly * p1, const gsl_poly * p2, gsl_poly * w)
{
    int i,j,k;
    int d1 = p1->degree;
    int d2 = p2->degree;
    double x;
    if (p1 -> size <= d1 + d2)
      {
	GSL_ERROR("size of destination poly not sufficient",GSL_EBADLEN);
      }
    if (w -> size  <= d1 + d2)
      {
	GSL_ERROR("size of workspace poly not sufficient",GSL_EBADLEN);
      }
    else
      {
	for (i = 0; i < w -> size; i++)
	  {
	    w->c[i] = 0.0;
	  }
	for (i = 0; i <= d1; i++)
	  {
	    x = p1->c[i];
	    for (j = 0; j <= d2; j++)
	      {
		k = i + j;
		w->c[k] += x * p2->c[j];
	      }
	  }
	for (i = 0; i <= d1+d2; i++)
	  {
	    p1->c[i] = w->c[i];
	  }
	p1->degree = d1 + d2;
      }
    return GSL_SUCCESS;
}

int gsl_poly_div(const gsl_poly * p1, const gsl_poly * p2, gsl_poly * q, gsl_poly * r)
{
    size_t d1 = p1->degree;
    size_t d2 = p2->degree;
    size_t dq = q->size - 1;
    size_t dr = r->size - 1;
    int k,j;

    if (d2 > d1) 
      {
	GSL_ERROR("degree of denominator should not exceed degree of numerator", GSL_EINVAL);
      }

    if (dq < d1)
      {
	GSL_ERROR("size of quotient not sufficient",GSL_EBADLEN);
      }

    if (dr < d1)
      {
	GSL_ERROR("size of remainder not sufficient",GSL_EBADLEN);
      }

    for (j = 0; j <= d1; j++) 
      {
	r->c[j] = p1->c[j];
	q->c[j] = 0.0;
      }

    for (k = d1 - d2; k >= 0; k--)
      {
	q->c[k] = r->c[d2+k]/p2->c[d2];
	for (j = d2 + k - 1; j >= k; j--)
	  {
	    r->c[j] -= q->c[k]*p2->c[j-k];
	  }
      }
    for (j = d2; j <= d1; j++)
      {
	r->c[j] = 0.0;
      } 
    gsl_poly_set_degree(q,GSL_DBL_EPSILON);
    gsl_poly_set_degree(r,GSL_DBL_EPSILON);
    return GSL_SUCCESS;
}

int gsl_poly_diff(const gsl_poly * p, gsl_poly * dp)
{
    int i;

    if (dp->size < p->degree)
      {
	GSL_ERROR("size of resulting poly not sufficient", GSL_EBADLEN);
      }

    for (i = 0; i < p->degree; i++)
      {
	dp->c[i] = (i+1) * p->c[i+1];
      }

    for (i = p->degree; i < dp->size; i++)
      {
	dp->c[i] = 0.0;
      }

    dp->degree = p->degree - 1;

    return GSL_SUCCESS;
}
