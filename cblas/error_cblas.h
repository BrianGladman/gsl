/* cblas/error_cblas.h
 *
 * Copyright (C) 2010 José Luis García Pallero
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

#ifndef __ERROR_CBLAS_H__
#define __ERROR_CBLAS_H__

#define CHECK_ARGS(FUNCTION,...) do { int pos = 0 ; \
    CBLAS_ERROR_##FUNCTION(pos,__VA_ARGS__) ; \
    if (pos) cblas_xerbla(pos,__FILE__,""); } while (0)

/* check if CBLAS_ORDER is correct */
#define CHECK_ORDER(pos,posIfError,order) \
if(((order)!=CblasRowMajor)&&((order)!=CblasColMajor)) \
     pos = posIfError;

/* check if CBLAS_TRANSPOSE is correct */
#define CHECK_TRANSPOSE(pos,posIfError,Trans) \
if(((Trans)!=CblasNoTrans)&&((Trans)!=CblasTrans)&&((Trans)!=CblasConjTrans)) \
    pos = posIfError;

/* check if CBLAS_UPLO is correct */
#define CHECK_UPLO(pos,posIfError,Uplo) \
if(((Uplo)!=CblasUpper)&&((Uplo)!=CblasLower)) \
    pos = posIfError;

/* check if CBLAS_DIAG is correct */
#define CHECK_DIAG(pos,posIfError,Diag) \
if(((Diag)!=CblasNonUnit)&&((Diag)!=CblasUnit)) \
    pos = posIfError;

/* check if CBLAS_SIDE is correct */
#define CHECK_SIDE(pos,posIfError,Side) \
if(((Side)!=CblasLeft)&&((Side)!=CblasRight)) \
    pos = posIfError;

/* check if a dimension argument is correct */
#define CHECK_DIM(pos,posIfError,dim) \
if((dim)<0) \
    pos = posIfError;

/* check if a stride argument is correct */
#define CHECK_STRIDE(pos,posIfError,stride) \
if((stride)==0) \
    pos = posIfError;

#endif /* __ERROR_CBLAS_H__ */
