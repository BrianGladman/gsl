/* err/env.c
 * 
 * Copyright (C) 2001, 2007 Brian Gough
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

#include <config.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>

#define STR_EQ(x,y) (strcmp((x),(y)) == 0)

void
gsl_err_env_setup (void)
{
  const char * p = getenv("GSL_ERROR_MODE") ;

  if (p == 0)  /* GSL_ERROR_MODE environment variable is not set */
    return ;

  if (*p == '\0') /* GSL_ERROR_MODE environment variable is empty */
    return ;

  printf("GSL_ERROR_MODE=\"") ;
  
  if (STR_EQ(p, "abort"))
    {
      gsl_set_error_handler (NULL);
      printf("abort") ;
    }

  printf("\"\n") ;
}





