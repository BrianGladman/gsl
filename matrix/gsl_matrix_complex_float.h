/* matrix/gsl_matrix_complex_float.h
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman, Brian Gough
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

#ifndef __GSL_MATRIX_COMPLEX_FLOAT_H__
#define __GSL_MATRIX_COMPLEX_FLOAT_H__

#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_vector_complex_float.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

typedef struct gsl_matrix_complex_float_struct gsl_matrix_complex_float;

struct gsl_matrix_complex_float_struct 
{
  size_t size1;
  size_t size2;
  size_t tda;
  float * data;
  gsl_block_complex_float * block;
} ;


gsl_matrix_complex_float * 
gsl_matrix_complex_float_alloc (const size_t n1, const size_t n2);

gsl_matrix_complex_float * 
gsl_matrix_complex_float_calloc (const size_t n1, const size_t n2);

gsl_matrix_complex_float * 
gsl_matrix_complex_float_alloc_from_block (gsl_block_complex_float * b, 
                                           const size_t offset, 
                                           const size_t n1, const size_t n2, const size_t d2);

gsl_matrix_complex_float * 
gsl_matrix_complex_float_alloc_from_matrix (gsl_matrix_complex_float * b,
                                            const size_t k1, const size_t k2,
                                            const size_t n1, const size_t n2);

gsl_vector_complex_float * 
gsl_vector_complex_float_alloc_row_from_matrix (gsl_matrix_complex_float * m,
                                                const size_t i);

gsl_vector_complex_float * 
gsl_vector_complex_float_alloc_col_from_matrix (gsl_matrix_complex_float * m,
                                                const size_t j);


void gsl_matrix_complex_float_free (gsl_matrix_complex_float * m);

gsl_matrix_complex_float gsl_matrix_complex_float_view (float * m, 
                                                        const size_t n1, 
                                                        const size_t n2);

void gsl_matrix_complex_float_set_zero (gsl_matrix_complex_float * m);
void gsl_matrix_complex_float_set_identity (gsl_matrix_complex_float * m);
void gsl_matrix_complex_float_set_all (gsl_matrix_complex_float * m, gsl_complex_float x);

gsl_complex_float * gsl_matrix_complex_float_ptr(const gsl_matrix_complex_float * m, const size_t i, const size_t j);
gsl_complex_float gsl_matrix_complex_float_get(const gsl_matrix_complex_float * m, const size_t i, const size_t j);
void gsl_matrix_complex_float_set(gsl_matrix_complex_float * m, const size_t i, const size_t j, const gsl_complex_float x);

int gsl_matrix_complex_float_fread (FILE * stream, gsl_matrix_complex_float * m) ;
int gsl_matrix_complex_float_fwrite (FILE * stream, const gsl_matrix_complex_float * m) ;
int gsl_matrix_complex_float_fscanf (FILE * stream, gsl_matrix_complex_float * m);
int gsl_matrix_complex_float_fprintf (FILE * stream, const gsl_matrix_complex_float * m, const char * format);

int gsl_matrix_complex_float_memcpy(gsl_matrix_complex_float * dest, const gsl_matrix_complex_float * src);
int gsl_matrix_complex_float_swap(gsl_matrix_complex_float * m1, const gsl_matrix_complex_float * m2);

int gsl_matrix_complex_float_swap_rows(gsl_matrix_complex_float * m, const size_t i, const size_t j);
int gsl_matrix_complex_float_swap_columns(gsl_matrix_complex_float * m, const size_t i, const size_t j);
int gsl_matrix_complex_float_swap_rowcol(gsl_matrix_complex_float * m, const size_t i, const size_t j);

int gsl_matrix_complex_float_transpose (gsl_matrix_complex_float * m);
int gsl_matrix_complex_float_transpose_memcpy (gsl_matrix_complex_float * dest, const gsl_matrix_complex_float * src);

gsl_matrix_complex_float gsl_matrix_complex_float_submatrix (gsl_matrix_complex_float * m, size_t i, size_t j, size_t n1, size_t n2);
gsl_vector_complex_float gsl_matrix_complex_float_row (gsl_matrix_complex_float * m, size_t i);
gsl_vector_complex_float gsl_matrix_complex_float_column (gsl_matrix_complex_float * m, size_t j);
gsl_vector_complex_float gsl_matrix_complex_float_diagonal (gsl_matrix_complex_float * m);

int gsl_matrix_complex_float_isnull (const gsl_matrix_complex_float * m);

/***********************************************************************/
/* The functions below are obsolete                                    */
/***********************************************************************/
int gsl_matrix_complex_float_get_row(gsl_vector_complex_float * v, const gsl_matrix_complex_float * m, const size_t i);
int gsl_matrix_complex_float_get_col(gsl_vector_complex_float * v, const gsl_matrix_complex_float * m, const size_t j);
int gsl_matrix_complex_float_set_row(gsl_matrix_complex_float * m, const size_t i, const gsl_vector_complex_float * v);
int gsl_matrix_complex_float_set_col(gsl_matrix_complex_float * m, const size_t j, const gsl_vector_complex_float * v);

int gsl_vector_complex_float_view_row_from_matrix (gsl_vector_complex_float * v, gsl_matrix_complex_float * m, const size_t i);
int gsl_vector_complex_float_view_col_from_matrix (gsl_vector_complex_float * v, gsl_matrix_complex_float * m, const size_t j);

int gsl_matrix_complex_float_view_from_vector (gsl_matrix_complex_float * m, gsl_vector_complex_float * base, const size_t offset, const size_t n1, const size_t n2, const size_t d2);

extern int gsl_check_range ;

#ifdef HAVE_INLINE

extern inline 
gsl_complex_float
gsl_matrix_complex_float_get(const gsl_matrix_complex_float * m, 
		     const size_t i, const size_t j)
{
  const gsl_complex_float zero = {{0,0}};

#ifndef GSL_RANGE_CHECK_OFF
  if (i >= m->size1)
    {
      GSL_ERROR_VAL("first index out of range", GSL_EINVAL, zero) ;
    }
  else if (j >= m->size2)
    {
      GSL_ERROR_VAL("second index out of range", GSL_EINVAL, zero) ;
    }
#endif
  return *(gsl_complex_float *)(m->data + 2*(i * m->tda + j)) ;
} 

extern inline 
void
gsl_matrix_complex_float_set(gsl_matrix_complex_float * m, 
		     const size_t i, const size_t j, const gsl_complex_float x)
{
#ifndef GSL_RANGE_CHECK_OFF
  if (i >= m->size1)
    {
      GSL_ERROR_VOID("first index out of range", GSL_EINVAL) ;
    }
  else if (j >= m->size2)
    {
      GSL_ERROR_VOID("second index out of range", GSL_EINVAL) ;
    }
#endif
  *(gsl_complex_float *)(m->data + 2*(i * m->tda + j)) = x ;
}
#endif /* HAVE_INLINE */

__END_DECLS

#endif /* __GSL_MATRIX_COMPLEX_FLOAT_H__ */
