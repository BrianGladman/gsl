/* blas/test_cases.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
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


#include "test_cases.h"

const float   c_2[] = { 2.0, 0.0 };
const float   c_3[] = { 3.0, 0.0 };
const double  z_2[] = { 2.0, 0.0 };
const double  z_3[] = { 3.0, 0.0 };

const float   vector_4_zero_f[] = { 0.0, 0.0, 0.0, 0.0 };
const double  vector_4_zero_d[] = { 0.0, 0.0, 0.0, 0.0 };

const float   vector_4_zero_c[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
const double  vector_4_zero_z[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

const float   vector_4_f[] = { -2.0, -1.0, 0.0, 3.0 };
const double  vector_4_d[] = { -2.0, -1.0, 0.0, 3.0 };

const float vector_4_c[] = 
{
  -2.0, 1.0, 
  -1.0, 1.0,
   0.0, 3.0,
   3.0, 5.0
};

const double vector_4_z[] = 
{
  -2.0, 1.0, 
  -1.0, 1.0,
   0.0, 3.0,
   3.0, 5.0
};

const float matrix_gen_4_f[] = 
{
  1.0,  0.0, -2.0, 1.0, 
  0.5,  3.0,  5.0, 1.0,
  2.0, -2.0, -1.0, 0.5,
  1.0, -1.0,  0.5, 0.0
};

const double matrix_gen_4_d[] = 
{
  1.0,  0.0, -2.0, 1.0, 
  0.5,  3.0,  5.0, 1.0,
  2.0, -2.0, -1.0, 0.5,
  1.0, -1.0,  0.5, 0.0
};

const float matrix_gen_4_c[] =
{
  1.0, 0.0,   0.0, 0.0,  -2.0, 1.0,   1.0, 1.0,
  0.5, 0.0,   3.0, 1.0,   5.0, 1.0,   0.0, 2.0,
  2.0, 1.0,  -2.0, 0.0,  -1.0, 0.0,   0.5, 3.0,
  1.0, 1.0,  -1.0, 2.0,   0.5, 1.0,   0.0, 1.0
};

const double matrix_gen_4_z[] =
{
  1.0, 0.0,   0.0, 0.0,  -2.0, 1.0,   1.0, 1.0,
  0.5, 0.0,   3.0, 1.0,   5.0, 1.0,   0.0, 2.0,
  2.0, 1.0,  -2.0, 0.0,  -1.0, 0.0,   0.5, 3.0,
  1.0, 1.0,  -1.0, 2.0,   0.5, 1.0,   0.0, 1.0
};

const float matrix_t_4_f[] = 
{
  1.0,  0.0, -2.0, 1.0, 
  0.5,  3.0,  5.0, 1.0,
  2.0, -2.0, -1.0, 0.5,
  1.0, -1.0,  0.5, 1.0
};

const double matrix_t_4_d[] = 
{
  1.0,  0.0, -2.0, 1.0, 
  0.5,  3.0,  5.0, 1.0,
  2.0, -2.0, -1.0, 0.5,
  1.0, -1.0,  0.5, 1.0
};

const float matrix_t_4_fpu[] = 
{
  1.0,  0.0, -2.0, 1.0, 
        3.0,  5.0, 1.0,
             -1.0, 0.5,
                   1.0
};

const float matrix_t_4_fpl[] = 
{
  1.0, 
  0.5,  3.0,
  2.0, -2.0, -1.0,
  1.0, -1.0,  0.5, 1.0
};

const double matrix_t_4_dpu[] = 
{
  1.0,  0.0, -2.0, 1.0, 
        3.0,  5.0, 1.0,
             -1.0, 0.5,
                   1.0
};

const double matrix_t_4_dpl[] = 
{
  1.0,
  0.5,  3.0,
  2.0, -2.0, -1.0,
  1.0, -1.0,  0.5, 1.0
};

const float matrix_gen_4_cpu[] =
{
  1.0, 0.0,   0.0, 0.0,  -2.0, 1.0,   1.0, 1.0,
              3.0, 1.0,   5.0, 1.0,   0.0, 2.0,
                         -1.0, 0.0,   0.5, 3.0,
                                      0.0, 1.0
};

const float matrix_gen_4_cpl[] =
{
  1.0, 0.0,
  0.5, 0.0,   3.0, 1.0,
  2.0, 1.0,  -2.0, 0.0,  -1.0, 0.0,
  1.0, 1.0,  -1.0, 2.0,   0.5, 1.0,   0.0, 1.0
};

const double matrix_gen_4_zpu[] =
{
  1.0, 0.0,   0.0, 0.0,  -2.0, 1.0,   1.0, 1.0,
              3.0, 1.0,   5.0, 1.0,   0.0, 2.0,
                         -1.0, 0.0,   0.5, 3.0,
                                      0.0, 1.0
};


const double matrix_gen_4_zpl[] =
{
  1.0, 0.0,
  0.5, 0.0,   3.0, 1.0,
  2.0, 1.0,  -2.0, 0.0,  -1.0, 0.0,
  1.0, 1.0,  -1.0, 2.0,   0.5, 1.0,   0.0, 1.0
};

const float matrix_her_4_c[] =
{
  1.0,  0.0,   0.0,  0.0,  -2.0,  1.0,   1.0, 1.0,
  0.0,  0.0,   3.0,  0.0,   5.0,  1.0,   0.0, 2.0,
 -2.0, -1.0,   5.0, -1.0,  -1.0,  0.0,   0.5, 3.0,
  1.0, -1.0,   0.0, -2.0,   0.5, -3.0,   1.0, 0.0
};

const double matrix_her_4_z[] =
{
  1.0,  0.0,   0.0,  0.0,  -2.0,  1.0,   1.0, 1.0,
  0.0,  0.0,   3.0,  0.0,   5.0,  1.0,   0.0, 2.0,
 -2.0, -1.0,   5.0, -1.0,  -1.0,  0.0,   0.5, 3.0,
  1.0, -1.0,   0.0, -2.0,   0.5, -3.0,   1.0, 0.0
};


const float matrix_her_4_cpu[] =
{
  1.0,  0.0,   0.0,  0.0,   -2.0,  1.0,   1.0, 1.0,
               3.0,  0.0,    5.0,  1.0,   0.0, 2.0,
                            -1.0,  0.0,   0.5, 3.0,
                                          1.0, 0.0
};

const double matrix_her_4_zpu[] =
{
  1.0,  0.0,   0.0,  0.0,   -2.0,  1.0,   1.0, 1.0,
               3.0,  0.0,    5.0,  1.0,   0.0, 2.0,
                            -1.0,  0.0,   0.5, 3.0,
                                          1.0, 0.0
};

const float matrix_her_4_cpl[] =
{
  1.0,  0.0, 
  0.0,  0.0,   3.0,   0.0, 
 -2.0, -1.0,   5.0,  -1.0,   -1.0,  0.0,
  1.0, -1.0,   0.0,  -2.0,    0.5, -3.0,   1.0, 0.0
};

const double matrix_her_4_zpl[] =
{
  1.0,  0.0, 
  0.0,  0.0,   3.0,   0.0, 
 -2.0, -1.0,   5.0,  -1.0,   -1.0,  0.0,
  1.0, -1.0,   0.0,  -2.0,    0.5, -3.0,   1.0, 0.0
};

const float   vector_7_f[] = { -3.0, 9.0, 2.5, 7.0, -11.0, 2.0, 5.0};
