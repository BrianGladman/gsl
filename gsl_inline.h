/* gsl_inline.h
 * 
 * Copyright (C) 2008 Brian Gough
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

#ifndef __GSL_INLINE_H__
#define __GSL_INLINE_H__

#ifdef HAVE_INLINE
#  ifdef __GNUC_STDC_INLINE__
#    define INLINE_DECL inline
#    define INLINE_FUN inline
#  else
#    define INLINE_DECL /* */
#    define INLINE_FUN extern inline
#  endif
#else
#  define INLINE_DECL /* */
#endif

#define GSL_RANGE_COND(x) (x)

#endif /* __GSL_INLINE_H__ */
