#ifndef GSL_ERRNO_H
#define GSL_ERRNO_H

#include <stdio.h>
#include <errno.h>

enum { 
  GSL_EDOM     = 1,   /* input domain error, e.g sqrt(-1) */
  GSL_ERANGE   = 2,   /* output range error, e.g. exp(1e100) */
  GSL_EFAULT   = 3,   /* invalid pointer */
  GSL_EINVAL   = 4,   /* invalid argument supplied by user */
  GSL_EFAILED  = 5,   /* generic failure */
  GSL_EFACTOR  = 6,   /* factorization failed */
  GSL_ESANITY  = 7,   /* sanity check failed - shouldn't happen */
  GSL_ENOMEM   = 8,   /* malloc failed */
  GSL_EBADFUNC = 9,   /* problem with user-supplied function */
  GSL_ERUNAWAY = 10,  /* iterative process is out of control */
  GSL_EMAXITER = 11,  /* exceeded max number of iterations */
  GSL_EZERODIV = 12,  /* tried to divide by zero */
  GSL_EBADTOL  = 13,  /* user specified an invalid tolerance */
  GSL_ETOL     = 14,  /* failed to reach the specified tolerance */
  GSL_EUNDRFLW = 15,  /* underflow */
  GSL_EOVRFLW  = 16,  /* overflow  */
  GSL_ELOSS    = 17,  /* loss of accuracy */
  GSL_EROUND   = 18,  /* failed because of roundoff error */
  GSL_EBADLEN  = 19,  /* matrix, vector lengths are not conformant */
  GSL_ESING    = 20,  /* apparent singularity detected */
  GSL_EUNSUP   = 21,  /* requested feature is not supported by the hardware */
  GSL_EUNIMPL  = 22   /* requested feature not (yet) implemented */
} ;

/* just to make things slightly clearer */

enum {
  GSL_SUCCESS = 0, 
  GSL_FAILURE = -1
};

void gsl_error (const char * reason, const char * file, int line,
		int gsl_errno);

void gsl_warning (const char * reason, const char * file, int line,
		  int gsl_errno) ;

void gsl_stream_printf (const char *label, const char *file,
			int line, const char *reason);

const char * gsl_strerror (int gsl_errno);

typedef void gsl_error_handler_t (const char * reason, const char * file,
				  int line, int gsl_errno);

typedef void gsl_stream_handler_t (const char * label, const char * file,
				   int line, const char * reason);

gsl_error_handler_t * 
gsl_set_error_handler (gsl_error_handler_t * new_handler);

gsl_stream_handler_t * 
gsl_set_stream_handler (gsl_stream_handler_t * new_handler);

FILE * gsl_set_stream (FILE * new_stream);

/* GSL_ERROR: call the error handler, and return the error code */

#define GSL_ERROR(reason, gsl_errno) \
       do { \
       gsl_error (reason, __FILE__, __LINE__, gsl_errno) ; \
       return gsl_errno ; \
       } while (0)

/* GSL_ERROR_RETURN: call the error handler, and return the given value */

#define GSL_ERROR_RETURN(reason, gsl_errno, value) \
       do { \
       gsl_error (reason, __FILE__, __LINE__, gsl_errno) ; \
       return value ; \
       } while (0)

/* GSL_ERROR_RETURN_NOTHING: call the error handler, and then return
   (for void functions which still need to generate an error) */

#define GSL_ERROR_RETURN_NOTHING(reason, gsl_errno) \
       do { \
       gsl_error (reason, __FILE__, __LINE__, gsl_errno) ; \
       return ; \
       } while (0)

/* GSL_ERROR_CONTINUE: call the error handler, and continue */

#define GSL_ERROR_CONTINUE(reason, gsl_errno) \
       do { \
       gsl_error (reason, __FILE__, __LINE__, gsl_errno) ; \
       } while (0)


/* GSL library code can occasionally generate warnings, which are not
   intended to be fatal. You can compile a version of the library with
   warnings turned off globally by defining the preprocessor constant
   GSL_WARNINGS_OFF. This turns off the warnings, but does not disable
   error handling in any way or turn off error messages.
 
   GSL_WARNING() is not intended for use in client code -- use
   GSL_MESSAGE() instead.  */
    
#ifdef GSL_WARNINGS_OFF   /* throw away warnings */
#define GSL_WARNING(warning, gsl_errno) \
       do { } while(0)
#else                     /* output all warnings */
#define GSL_WARNING(warning, gsl_errno) \
       do { \
       gsl_warning (warning, __FILE__, __LINE__, gsl_errno) ; \
       } while (0)
#endif

/* Warnings can also be turned off at runtime by setting the variable
   gsl_warnings_off to a non-zero value */

extern int gsl_warnings_off ;

#endif /* GSL_ERRNO_H */

