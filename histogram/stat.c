/* gsl_histogram_stat.c
 * Copyright (C) 2000  Simone Piccardi
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */
/***************************************************************
 *
 * File gsl_histogram_stat.c: 
 * Routines for statisticalcomputations on histograms. 
 * Need GSL library and header.
 * Contains the routines:
 * gsl_histogram_mean    compute histogram mean
 * gsl_histogram_sigma   compute histogram sigma
 *
 * Author: S. Piccardi
 * Jan. 2000
 *
 * $Id$
 *
 ***************************************************************/
#include <config.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_histogram.h>

/* FIXME: these algorithms could be made more stable by using the
   recurrences found in the statistics directory. Also we should
   probably skip negative values in the histogram h->bin[i] < 0, since
   those correspond to negative weights (BJG) */

double
gsl_histogram_mean (const gsl_histogram * h)
{
  const size_t n = h->n;
  size_t i;
  double meansum = 0;
  double entries = 0;

  for (i = 0; i < n; i++)
    {
      double binvalue = (h->range[i + 1] - h->range[i]) / 2;
      meansum += binvalue * h->bin[i];
      entries += h->bin[i];
    }

  return meansum / entries;
}

double
gsl_histogram_sigma (const gsl_histogram * h)
{
  const size_t n = h->n;
  size_t i;
  double meansum = 0;
  double sigmasum = 0;
  double entries = 0;

  for (i = 0; i < n; i++)
    {
      double binvalue = ((h->range[i + 1]) - (h->range[i])) / 2;
      sigmasum += binvalue * binvalue * h->bin[i];
      meansum += binvalue * h->bin[i];
      entries += h->bin[i];
    }

  return sigmasum / entries - meansum / entries * meansum / entries;
}
