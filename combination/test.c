/* combination/test.c
 * based on permutation/test.c by Brian Gough
 * 
 * Copyright (C) 2001 Szymon Jaroszewicz
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
#include <gsl/gsl_combination.h>
#include <gsl/gsl_test.h>

size_t c63[20][3] = {
  { 0, 1, 2 },  { 0, 1, 3 },  { 0, 1, 4 },  { 0, 1, 5 },
  { 0, 2, 3 },  { 0, 2, 4 },  { 0, 2, 5 },  { 0, 3, 4 },
  { 0, 3, 5 },  { 0, 4, 5 },  { 1, 2, 3 },  { 1, 2, 4 },
  { 1, 2, 5 },  { 1, 3, 4 },  { 1, 3, 5 },  { 1, 4, 5 },
  { 2, 3, 4 },  { 2, 3, 5 },  { 2, 4, 5 },  { 3, 4, 5 }
} ;


int 
main (void)
{
  int i = 0, j, status = 0;
  gsl_combination * c ;

  c = gsl_combination_alloc (6,3);

  gsl_combination_init_first (c);
  
  do 
    {
      if ( i >= 20 )
        {
	  status = 1;
          break;
	}
      for (j = 0; j < 3; j++)
        {
          status |= (c->data[j] != c63[i][j]);
        }
      i++;
    }
  while (gsl_combination_next(c) == GSL_SUCCESS);

  gsl_test(status, "gsl_combination_next, 6 choose 3 combination, 20 steps");

  gsl_combination_next(c);
  gsl_combination_next(c);
  gsl_combination_next(c);
  for (j = 0; j < 3; j++)
    {
      status |= (c->data[j] != c63[19][j]);
    }
  gsl_test(status, "gsl_combination_next on the last combination");


  gsl_combination_init_last (c);

  i = 19;
  do 
    {
      if ( i < 0 )
        {
	  status = 1;
          break;
	}
      for (j = 0; j < 3; j++)
        {
          status |= (c->data[j] != c63[i][j]);
        }
      i--;
    }
  while (gsl_combination_prev(c) == GSL_SUCCESS);

  gsl_test(status, "gsl_combination_prev, 6 choose 3 combination, 20 steps");

  gsl_combination_prev(c);
  gsl_combination_prev(c);
  gsl_combination_prev(c);
  for (j = 0; j < 3; j++)
    {
      status |= (c->data[j] != c63[0][j]);
    }
  gsl_test(status, "gsl_combination_prev on the first combination");
  gsl_combination_free (c);

  c = gsl_combination_calloc(7, 0);
  /* should return GSL_FAILURE every time */
  status |= (gsl_combination_next(c) != GSL_FAILURE);
  status |= (gsl_combination_next(c) != GSL_FAILURE);
  status |= (gsl_combination_prev(c) != GSL_FAILURE);
  status |= (gsl_combination_prev(c) != GSL_FAILURE);
  gsl_test(status, "gsl_combination 7 choose 0");
  gsl_combination_free (c);

  c = gsl_combination_calloc(7, 7);
  /* should return GSL_FAILURE every time */
  for(j = 0; j < 7; j++)
  {
    status |= (gsl_combination_get(c, j) != j);
  }
  status |= (gsl_combination_next(c) != GSL_FAILURE);
  for(j = 0; j < 7; j++)
  {
    status |= (gsl_combination_get(c, j) != j);
  }
  status |= (gsl_combination_next(c) != GSL_FAILURE);
  for(j = 0; j < 7; j++)
  {
    status |= (gsl_combination_get(c, j) != j);
  }
  status |= (gsl_combination_prev(c) != GSL_FAILURE);
  for(j = 0; j < 7; j++)
  {
    status |= (gsl_combination_get(c, j) != j);
  }
  status |= (gsl_combination_prev(c) != GSL_FAILURE);
  for(j = 0; j < 7; j++)
  {
    status |= (gsl_combination_get(c, j) != j);
  }
  gsl_test(status, "gsl_combination 7 choose 7");
  gsl_combination_free (c);

  exit (gsl_test_summary());
}
