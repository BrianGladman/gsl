/* permutation/test.c
 * 
 * Copyright (C) 2000 Brian Gough
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
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_permute_double.h>
#include <gsl/gsl_test.h>

unsigned int p5[120][5] = {
  {0, 1, 2, 3, 4}, {0, 1, 2, 4, 3}, {0, 1, 3, 2, 4}, {0, 1, 3, 4, 2},
  {0, 1, 4, 2, 3}, {0, 1, 4, 3, 2}, {0, 2, 1, 3, 4}, {0, 2, 1, 4, 3},
  {0, 2, 3, 1, 4}, {0, 2, 3, 4, 1}, {0, 2, 4, 1, 3}, {0, 2, 4, 3, 1},
  {0, 3, 1, 2, 4}, {0, 3, 1, 4, 2}, {0, 3, 2, 1, 4}, {0, 3, 2, 4, 1},
  {0, 3, 4, 1, 2}, {0, 3, 4, 2, 1}, {0, 4, 1, 2, 3}, {0, 4, 1, 3, 2},
  {0, 4, 2, 1, 3}, {0, 4, 2, 3, 1}, {0, 4, 3, 1, 2}, {0, 4, 3, 2, 1},
  {1, 0, 2, 3, 4}, {1, 0, 2, 4, 3}, {1, 0, 3, 2, 4}, {1, 0, 3, 4, 2},
  {1, 0, 4, 2, 3}, {1, 0, 4, 3, 2}, {1, 2, 0, 3, 4}, {1, 2, 0, 4, 3},
  {1, 2, 3, 0, 4}, {1, 2, 3, 4, 0}, {1, 2, 4, 0, 3}, {1, 2, 4, 3, 0},
  {1, 3, 0, 2, 4}, {1, 3, 0, 4, 2}, {1, 3, 2, 0, 4}, {1, 3, 2, 4, 0},
  {1, 3, 4, 0, 2}, {1, 3, 4, 2, 0}, {1, 4, 0, 2, 3}, {1, 4, 0, 3, 2},
  {1, 4, 2, 0, 3}, {1, 4, 2, 3, 0}, {1, 4, 3, 0, 2}, {1, 4, 3, 2, 0},
  {2, 0, 1, 3, 4}, {2, 0, 1, 4, 3}, {2, 0, 3, 1, 4}, {2, 0, 3, 4, 1},
  {2, 0, 4, 1, 3}, {2, 0, 4, 3, 1}, {2, 1, 0, 3, 4}, {2, 1, 0, 4, 3},
  {2, 1, 3, 0, 4}, {2, 1, 3, 4, 0}, {2, 1, 4, 0, 3}, {2, 1, 4, 3, 0},
  {2, 3, 0, 1, 4}, {2, 3, 0, 4, 1}, {2, 3, 1, 0, 4}, {2, 3, 1, 4, 0},
  {2, 3, 4, 0, 1}, {2, 3, 4, 1, 0}, {2, 4, 0, 1, 3}, {2, 4, 0, 3, 1},
  {2, 4, 1, 0, 3}, {2, 4, 1, 3, 0}, {2, 4, 3, 0, 1}, {2, 4, 3, 1, 0},
  {3, 0, 1, 2, 4}, {3, 0, 1, 4, 2}, {3, 0, 2, 1, 4}, {3, 0, 2, 4, 1},
  {3, 0, 4, 1, 2}, {3, 0, 4, 2, 1}, {3, 1, 0, 2, 4}, {3, 1, 0, 4, 2},
  {3, 1, 2, 0, 4}, {3, 1, 2, 4, 0}, {3, 1, 4, 0, 2}, {3, 1, 4, 2, 0},
  {3, 2, 0, 1, 4}, {3, 2, 0, 4, 1}, {3, 2, 1, 0, 4}, {3, 2, 1, 4, 0},
  {3, 2, 4, 0, 1}, {3, 2, 4, 1, 0}, {3, 4, 0, 1, 2}, {3, 4, 0, 2, 1},
  {3, 4, 1, 0, 2}, {3, 4, 1, 2, 0}, {3, 4, 2, 0, 1}, {3, 4, 2, 1, 0},
  {4, 0, 1, 2, 3}, {4, 0, 1, 3, 2}, {4, 0, 2, 1, 3}, {4, 0, 2, 3, 1},
  {4, 0, 3, 1, 2}, {4, 0, 3, 2, 1}, {4, 1, 0, 2, 3}, {4, 1, 0, 3, 2},
  {4, 1, 2, 0, 3}, {4, 1, 2, 3, 0}, {4, 1, 3, 0, 2}, {4, 1, 3, 2, 0},
  {4, 2, 0, 1, 3}, {4, 2, 0, 3, 1}, {4, 2, 1, 0, 3}, {4, 2, 1, 3, 0},
  {4, 2, 3, 0, 1}, {4, 2, 3, 1, 0}, {4, 3, 0, 1, 2}, {4, 3, 0, 2, 1},
  {4, 3, 1, 0, 2}, {4, 3, 1, 2, 0}, {4, 3, 2, 0, 1}, {4, 3, 2, 1, 0}
} ;

int 
main (void)
{
  int i = 0, j, status = 0;
  gsl_permutation * p ;

  p = gsl_permutation_alloc (5);

  gsl_permutation_init (p);
  
  do 
    {
      for (j = 0; j < 5; j++)
        {
          status |= (p->data[j] != p5[i][j]);
        }

      i++;
    }
  while (gsl_permutation_next(p) == GSL_SUCCESS);

  gsl_test(status, "gsl_permutation_next, 5-th order permutation, 120 steps");


  do 
    {
      i--;
      
      for (j = 0; j < 5; j++)
        {
          status |= (p->data[j] != p5[i][j]);
        }
    }
  while (gsl_permutation_prev(p) == GSL_SUCCESS);

  gsl_test(status, "gsl_permutation_prev, 5-th order permutation, 120 steps");

  gsl_permutation_free (p);

  {
    int i;
    int status = 0 ;

    gsl_permutation * p1 = gsl_permutation_alloc (5);
    gsl_permutation * p2 = gsl_permutation_alloc (5);
    gsl_permutation * p = gsl_permutation_alloc (5);

    double v[5] = { 100.0, 101.0, 102.0, 103.0, 104.0 } ;

    gsl_permutation_init (p1);

    do 
      {
        gsl_permutation_init (p2);

        do 
          {
            double x[5], y[5];

            /* Compute x= p1 p2 v */
            memcpy (x, v, 5 * sizeof(double));
            gsl_permute (p2->data, x, 1, 5);
            gsl_permute (p1->data, x, 1, 5);

            /* Compute P= p1 p2, y = P v */
            gsl_permutation_mul (p, p1, p2);
            memcpy (y, v, 5 * sizeof(double));
            gsl_permute (p->data, y, 1, 5);

            for (i = 0; i < 5; i++) 
              {
                if (x[i] != y[i]) 
                  status = 1;
              }

            if (status == 1)
              break;
          }
        while (gsl_permutation_next(p2) == GSL_SUCCESS);
        
        if (status == 1)
          break;
      }
    while (gsl_permutation_next(p1) == GSL_SUCCESS);
    
    gsl_permutation_free (p1);
    gsl_permutation_free (p2);
    gsl_permutation_free (p);

    gsl_test(status, "gsl_permutation_mul, all 5-th order combinations");
  }

  exit (gsl_test_summary());
}
