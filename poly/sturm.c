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



gsl_poly_sturm *
gsl_poly_sturm_alloc (const size_t size)
{
  int i = 0, k = 0;
  gsl_poly_sturm *ss = (gsl_poly_sturm *) malloc (sizeof (gsl_poly_sturm));

  if (ss == 0)
    {
      GSL_ERROR_VAL ("failed to allocate gsl_poly_sturm struct", GSL_ENOMEM,
		     0);
    }

  ss->sturmseq = malloc (size * sizeof (gsl_poly *));

  if (ss->sturmseq == 0)
    {
      free (ss);
      GSL_ERROR_VAL ("failed to allocate gsl_poly_sturm data array",
		     GSL_ENOMEM, 0);
    }

  for (i = 0; i < size; i++)
    {
      ss->sturmseq[i] = gsl_poly_alloc (size);

      if (ss->sturmseq[i] == 0)
	{
	  for (k = 0; k < i; k++)
	    {
	      gsl_poly_free (ss->sturmseq[k]);
	    }
	  free (ss);
	  GSL_ERROR_VAL ("failed to allocate gsl_poly_sturm data array",
			 GSL_ENOMEM, 0);

	}
    }
  ss->size = size;

  return ss;
}

gsl_poly_sturm *
gsl_poly_sturm_calloc (const size_t size)
{
  int i = 0, k = 0;
  gsl_poly_sturm *ss = (gsl_poly_sturm *) malloc (sizeof (gsl_poly_sturm));

  if (ss == 0)
    {
      GSL_ERROR_VAL ("failed to allocate gsl_poly_sturm struct", GSL_ENOMEM,
		     0);
    }

  ss->sturmseq = malloc (size * sizeof (gsl_poly *));

  if (ss->sturmseq == 0)
    {
      free (ss);
      GSL_ERROR_VAL ("failed to allocate gsl_poly_sturm data array",
		     GSL_ENOMEM, 0);
    }

  for (i = 0; i < size; i++)
    {
      ss->sturmseq[i] = gsl_poly_calloc (size);

      if (ss->sturmseq[i] == 0)
	{
	  for (k = 0; k < i; k++)
	    {
	      gsl_poly_free (ss->sturmseq[k]);
	    }
	  free (ss);
	  GSL_ERROR_VAL ("failed to allocate gsl_poly_sturm data array",
			 GSL_ENOMEM, 0);
	}
    }
  ss->size = size;

  return ss;
}

void
gsl_poly_sturm_free (gsl_poly_sturm * ss)
{
  int i = 0;
  int size = ss->size;

  for (i = 0; i < size; i++)
    {
      gsl_poly_free (ss->sturmseq[i]);
    }

  free (ss->sturmseq);

  free (ss);

}

int
gsl_poly_sturm_build (gsl_poly_sturm * ss, const gsl_poly * p, gsl_poly * w)
{
  int d = p->degree;
  int index = 0;
  double alpha = 0.0;

  gsl_poly *q = w;

  if (q->size <= p->degree)
    {
      GSL_ERROR ("Size of workspace poly not sufficient", GSL_EBADLEN);
    }


  alpha = gsl_poly_get (p, d);

  if (ss->size < p->degree)
    {
      free (q);
      GSL_ERROR ("Sturm sequence not large enough", GSL_EINVAL);
    }

  gsl_poly_memcpy (ss->sturmseq[0], p);

  gsl_poly_scale (ss->sturmseq[0], alpha);

  gsl_poly_diff (ss->sturmseq[0], ss->sturmseq[1]);

  alpha = fabs (gsl_poly_get (ss->sturmseq[1], d - 1));
  gsl_poly_scale (ss->sturmseq[1], 1.0 / alpha);

  index = 2;

  while (ss->sturmseq[index - 1]->degree > 0)
    {
      gsl_poly_div (ss->sturmseq[index - 2], ss->sturmseq[index - 1], q,
		    ss->sturmseq[index]);
      d = ss->sturmseq[index]->degree;
      alpha = fabs (gsl_poly_get (ss->sturmseq[index], d));
      if (ss->sturmseq[index]->degree > 0)
	{
	  gsl_poly_scale (ss->sturmseq[index], -1.0 / alpha);
	}
      else
	{
	  gsl_poly_scale (ss->sturmseq[index], -1.0);
	}
      index++;
    }

  ss->index = index - 1;
  return GSL_SUCCESS;
}

int
gsl_poly_sturm_changes (const gsl_poly_sturm * ss, double a)
{
  int changes = 0;
  double u = 0.0, v = 0.0;
  int i = 0;
  int index = ss->index;

  v = gsl_poly_eval2 (ss->sturmseq[0], a);

  for (i = 1; i <= index; i++)
    {
      u = gsl_poly_eval2 (ss->sturmseq[i], a);
      if (v == 0.0 || v * u < 0.0)
	{
	  changes++;
	}
      v = u;
    }


  return changes;

}

int
gsl_poly_sturm_numroots (const gsl_poly_sturm * ss, double a, double b)
{
  int nr = gsl_poly_sturm_changes (ss, a) - gsl_poly_sturm_changes (ss, b);

  return nr;
}

void
gsl_poly_sturm_dump (const gsl_poly_sturm * ss)
{

  int i = 0, j = 0;
  gsl_poly *p;

  for (i = 0; i <= ss->index; i++)
    {
      p = ss->sturmseq[i];
      for (j = p->degree; j >= 0; j--)
	{
	  printf ("%f ", p->c[j]);
	}
      printf ("\n");
    }
}

  
#ifdef STURM
  {
    gsl_poly *p = gsl_poly_calloc (6);
    int c1, c2, c3, nr;
    gsl_poly_sturm *pss = gsl_poly_sturm_calloc (6);
    gsl_poly *w = gsl_poly_calloc (6);
    gsl_poly_set (p, 0, -1);
    gsl_poly_set (p, 1, -3);
    gsl_poly_set (p, 5, 1);
    gsl_test (pss == 0, "gsl_poly_sturm_calloc returns valid pointer");
    gsl_poly_sturm_build (pss, p, w);
    gsl_test (pss->sturmseq[2]->degree != 1,
              "Sturm function 2, degree = 1");
    gsl_test (pss->index != 3, "Sturm function index = 3");
    gsl_test_rel (pss->sturmseq[2]->c[0], 10.0 / 24.0, 1.0e-9,
                  "Sturm function 2, x^0");
    gsl_test_rel (pss->sturmseq[2]->c[1], 1.0, 1.0e-9,
                  "Sturm function 2, x^1");
    gsl_test_rel (pss->sturmseq[3]->c[0], 5.69859182098765e-01, 1.0e-9,
                  "Sturm function 3, x^0");
    gsl_test (pss->sturmseq[3]->degree != 0,
              "Sturm function 3, degree = 0");
    gsl_test (pss->sturmseq[4]->degree != 0,
              "Sturm function 4, degree = 0");
    gsl_test (pss->sturmseq[5]->degree != 0,
              "Sturm function 5, degree = 0");
    c1 = gsl_poly_sturm_changes (pss, -2.0);
    c2 = gsl_poly_sturm_changes (pss, 0.0);
    c3 = gsl_poly_sturm_changes (pss, +2.0);
    gsl_test_rel (c1, 3.0, 1.0e-9, "Sturm changes in x = -2.0");
    gsl_test_rel (c2, 1.0, 1.0e-9, "Sturm changes in x =  0.0");
    gsl_test_abs (c3, 0.0, 1.0e-9, "Sturm changes in x = +2.0");
    nr = gsl_poly_sturm_numroots (pss, -2.0, 2.0);
    gsl_test_rel (nr, 3, 1.0e-9, "Number of roots in [-2,2]");
    
    free (p);
    free (w);
    free (pss);
  }
  
  {
    gsl_poly_sturm *pss = gsl_poly_sturm_calloc (11);
    gsl_poly *p = gsl_poly_calloc (11);
    gsl_poly *w = gsl_poly_calloc (11);
    double a = 0.0, b = 0.0;
    int i = 0, nr = 0;
    
    /* Wilkinson poly */
    gsl_poly_set (p, 0, 3628800.0);
    gsl_poly_set (p, 1, -10628640.0);
    gsl_poly_set (p, 2, 12753576.0);
    gsl_poly_set (p, 3, -8409500.0);
    gsl_poly_set (p, 4, 3416930.0);
    gsl_poly_set (p, 5, -902055.0);
    gsl_poly_set (p, 6, 157773.0);
    gsl_poly_set (p, 7, -18150.0);
    gsl_poly_set (p, 8, 1320.0);
    gsl_poly_set (p, 9, -55.0);
    gsl_poly_set (p, 10, 1.0);
    gsl_poly_sturm_build (pss, p, w);
    for (i = 0; i < 11; i++)
      {
        b = ((double) i) + 1.0 / 2.0;
        nr = gsl_poly_sturm_numroots (pss, a, b);
        gsl_test_rel (nr, i, 1.0e-12,
                      "Number of roots in [0,%f] Wilkinson p10");
      }
    
    gsl_poly_sturm_free (pss);
    gsl_poly_free (p);
    gsl_poly_free (w);
  }
  
  {
    gsl_poly_sturm *pss = gsl_poly_sturm_alloc (3);
    gsl_poly *p = gsl_poly_calloc (3);
    gsl_poly *w = gsl_poly_calloc (3);
    int nr = 0;
    
    gsl_poly_set (p, 0, 1.0);
    gsl_poly_set (p, 2, 0.825);
    
    gsl_poly_sturm_build (pss, p, w);
    nr = gsl_poly_sturm_numroots (pss, 0.0, 1.0);
    gsl_test_abs (nr, 0, 1.0e-12, "Number of roots of 1+0.825x^2 in [0,1]");
    
    gsl_poly_sturm_free (pss);
    gsl_poly_free (p);
    gsl_poly_free (w);
  }
#endif
