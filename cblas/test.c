/* blas/test.c
 * 
 * Copyright (C) 2001 Brian Gough
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

#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

void test_dot ();
void test_nrm2 ();
void test_asum ();
void test_amax ();
void test_axpy ();
void test_copy ();
void test_swap ();
void test_scal ();
void test_rotg ();
void test_rot ();
void test_rotmg ();
void test_rotm ();
void test_gemv ();
void test_gbmv ();
void test_trmv ();
void test_tbmv ();
void test_tpmv ();
void test_symv ();
void test_hemv ();
void test_hbmv ();
void test_sbmv ();
void test_hpmv ();
void test_spmv ();
void test_trsv ();
void test_tbsv ();
void test_tpsv ();
void test_ger ();
void test_syr ();
void test_her ();
void test_hpr ();
void test_spr ();
void test_syr2 ();
void test_spr2 ();
void test_her2 ();
void test_hpr2 ();
void test_gemm ();
void test_symm ();
void test_hemm ();
void test_syrk ();
void test_herk ();
void test_syr2k ();
void test_her2k ();
void test_trmm ();
void test_trsm ();

int 
main ()
{
  gsl_ieee_env_setup ();

#include "tests.c"

  return gsl_test_summary();
}

