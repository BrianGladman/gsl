#ifndef _GSL_ERRNO_H
#define _GSL_ERRNO_H

#include <errno.h>

#define GSL_EDOM       1   /* input domain error, e.g sqrt(-1) */
#define GSL_ERANGE     2   /* output range error, e.g. exp(1e100) */
#define GSL_EFAULT     3 
#define GSL_EINVAL     4 
#define GSL_EFAILED    5   /* generic failure */
#define GSL_EFACTOR    6   /* factorization failed */
#define GSL_ESANITY    7   /* check failed - shouldn't happen */
#define GSL_ENOMEM     8   /* malloc failed */
#define GSL_EBADFUNC   9   /* problem with user-supplied function */
#define GSL_ERUNAWAY   10  /* iterative process is out of control */
#define GSL_ETIMEOUT   11  /* exceeded max number of iterations */
#define GSL_EZERODIV   12  /* tried to divide by zero */
#define GSL_ETOL       13  /* user specified an invalid tolerance */
#define GSL_EUNDRFLW   14  /* underflow */
#define GSL_EOVRFLW    15  /* overflow  */
#define GSL_ELOSS      16  /* loss of accuracy */

/* just to make things slightly clearer */
enum {GSL_SUCCESS = 0, GSL_FAILURE = -1};

void gsl_error (const char *reason, const char *file, int line, int gsl_errno);
void gsl_warning(const char * reason, const char * file, int line) ;

typedef void gsl_errhandler_t (const char *reason,
			       const char *file,
			       int line,
			       int gsl_errno);

void gsl_empty_error_handler (const char *reason,
			      const char * file,
			      int line,
			      int gsl_errno);

gsl_errhandler_t * gsl_set_error_handler (gsl_errhandler_t * new_handler);

#ifdef GSL_THREAD_SAFE
#define GSL_ERRHANDLER_OFF
#endif

/* call error handler, returning error status in POSIX fashion */

#ifdef GSL_ERRHANDLER_OFF
#define GSL_ERROR(reason, gsl_errno) \
       do { \
       return gsl_errno ; \
       } while (0)
/* call error handler, returning a specified value [non-POSIX behaviour] */
#define GSL_ERROR_RETURN(reason, gsl_errno, value) \
       do { \
       GSL_WARNING(reason) ; \
       return value ; \
       } while (0)
#else
#define GSL_ERROR(reason, gsl_errno) \
       do { \
       gsl_error (reason, __FILE__, __LINE__, gsl_errno) ; \
       return gsl_errno ; \
       } while (0)
/* call error handler, returning a specified value [non-POSIX behaviour] */
#define GSL_ERROR_RETURN(reason, gsl_errno, value) \
       do { \
       GSL_WARNING(reason) ; \
       gsl_error (reason, __FILE__, __LINE__, gsl_errno) ; \
       return value ; \
       } while (0)
#endif /* GSL_ERRHANDLER_OFF */

/* GSL library code can occasionally generate warnings, which are not
   intended to be fatal. You can compile a version of the library with
   warnings turned off globally by defining the preprocessor constant
   GSL_WARNINGS_OFF. This turns off the warnings, but does not disable
   error handling in any way or turn off error messages.
 
   GSL_WARNING() is not intended for use in client code -- use
   GSL_MESSAGE() instead.  */

extern int gsl_warnings_off ;

/* Warnings can also be turned off at runtime by setting the variable
   gsl_warnings_off to a non-zero value */
     
#ifdef GSL_WARNINGS_OFF   /* throw away warnings */
#define GSL_WARNING(warning) do { } while(0)
#else                     /* output all warnings */
#define GSL_WARNING(warning) \
       do { \
       gsl_warning (warning, __FILE__, __LINE__) ; \
       } while (0)
#endif

#endif /* !_GSL_ERRNO_H */
