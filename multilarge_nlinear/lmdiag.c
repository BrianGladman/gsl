/* multilarge_nlinear/lmdiag.c
 * 
 * Copyright (C) 2015 Patrick Alken
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/*
 * This module handles the updating of the scaling matrix D
 * in the solution step:
 *
 * (J^T J + mu D^T D) dx = -J^T f
 *
 * according to several different strategies.
 */

static int init_diag_levenberg(const gsl_matrix * JTJ, gsl_vector * DTD);
static int update_diag_levenberg(const gsl_matrix * JTJ,
                                 gsl_vector * DTD);

static int init_diag_marquardt(const gsl_matrix * JTJ, gsl_vector * DTD);
static int update_diag_marquardt (const gsl_matrix * JTJ,
                                  gsl_vector * DTD);

static int init_diag_more(const gsl_matrix * JTJ, gsl_vector * DTD);
static int update_diag_more(const gsl_matrix * JTJ, gsl_vector * DTD);

/* Levenberg scaling, D = I */
static int
init_diag_levenberg(const gsl_matrix * JTJ, gsl_vector * DTD)
{
  (void)JTJ; /* avoid unused parameter warning */
  gsl_vector_set_all(DTD, 1.0);
  return GSL_SUCCESS;
}

static int
update_diag_levenberg(const gsl_matrix * JTJ, gsl_vector * DTD)
{
  (void)JTJ; /* avoid unused parameter warning */
  (void)DTD; /* avoid unused parameter warning */

  /* nothing to do */
  return GSL_SUCCESS;
}

/* initialize diagonal scaling matrix D according to Marquardt method */
static int
init_diag_marquardt(const gsl_matrix * JTJ, gsl_vector * DTD)
{
  return update_diag_marquardt(JTJ, DTD);
}

/* update diagonal scaling matrix D according to Marquardt method */
static int
update_diag_marquardt (const gsl_matrix * JTJ, gsl_vector * DTD)
{
  const size_t p = JTJ->size2;
  size_t j;

  for (j = 0; j < p; j++)
    {
      double ajj = gsl_matrix_get(JTJ, j, j);

      if (ajj == 0.0)
        ajj = 1.0;

      gsl_vector_set(DTD, j, ajj);
    }

  return GSL_SUCCESS;
}

/* initialize diagonal scaling matrix D according to Eq 6.3 of
 * More, 1978 */
static int
init_diag_more(const gsl_matrix * JTJ, gsl_vector * DTD)
{
  int status;

  gsl_vector_set_zero(DTD);
  status = update_diag_more(JTJ, DTD);

  return status;
}

/* update diagonal scaling matrix D according to Eq. 6.3 of
 * More, 1978 */
static int
update_diag_more (const gsl_matrix * JTJ, gsl_vector * DTD)
{
  const size_t p = JTJ->size2;
  size_t j;

  for (j = 0; j < p; j++)
    {
      double ajj = gsl_matrix_get(JTJ, j, j);
      double *dj = gsl_vector_ptr(DTD, j);

      if (ajj == 0.0)
        ajj = 1.0;

      *dj = GSL_MAX(*dj, ajj);
    }

  return GSL_SUCCESS;
}

static const gsl_multilarge_nlinear_scale levenberg_type =
{
  "levenberg",
  init_diag_levenberg,
  update_diag_levenberg
};

static const gsl_multilarge_nlinear_scale marquardt_type =
{
  "marquardt",
  init_diag_marquardt,
  update_diag_marquardt
};

static const gsl_multilarge_nlinear_scale more_type =
{
  "more",
  init_diag_more,
  update_diag_more
};

const gsl_multilarge_nlinear_scale *gsl_multilarge_nlinear_scale_levenberg = &levenberg_type;
const gsl_multilarge_nlinear_scale *gsl_multilarge_nlinear_scale_marquardt = &marquardt_type;
const gsl_multilarge_nlinear_scale *gsl_multilarge_nlinear_scale_more = &more_type;
