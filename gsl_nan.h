/* gsl_nan.h
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

#ifndef __GSL_NAN_H__
#define __GSL_NAN_H__

#ifdef NAN
#define GSL_NAN NAN
#define GSL_POSINF HUGE_VAL
#define GSL_NEGINF (-HUGE_VAL)
#define GSL_POSZERO (+0)
#define GSL_NEGZERO (-0)
#elif defined(_MSC_VER) /* Microsoft Visual C++ */
#define GSL_NAN _FPCLASS_QNAN
#define GSL_POSINF _FPCLASS_PINF
#define GSL_NEGINF _FPCLASS_NINF
#define GSL_POSZERO _FPCLASS_PZ
#define GSL_NEGZERO _FPCLASS_NZ
#else
#define GSL_NAN (0.0/0.0)
#define GSL_POSINF (+1.0/0.0)
#define GSL_NEGINF (-1.0/0.0)
#define GSL_POSZERO (+0)
#define GSL_NEGZERO (-0)
#endif

#endif /* __GSL_NAN_H__ */
