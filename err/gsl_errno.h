#ifndef _GSL_ERRNO_H
#define _GSL_ERRNO_H

#include <errno.h>

#define GSL_EDOM       1  /* domain error for input, e.g. asin(10),sqrt(-1) */
#define GSL_ERANGE     2  /* range error for output, e.g. exp(1e100) */
#define GSL_EFAULT     3
#define GSL_EINVAL     4
#define GSL_EFAILED    5  /* generic failure */
#define GSL_EFACTOR    6  /* factorization failed */
#define GSL_ESANITY    7  /* sanity check failed -- this shouldn't happen */
#define GSL_ENOMEM     8  /* malloc failed */
#define GSL_EBADFUNC   9  /* called a user-supplied function and it didn't
                             work properly, e.g. we needed a real number and
                             it gave us NAN */
#define GSL_ERUNAWAY  10  /* the next iteration was going to do something
                             ridiculous, e.g. a root finder was close to an
                             extremum and about to make a huge jump */
#define GSL_ETIMEOUT  11  /* exceeded the allowed number of iterations or took
                             too much time doing something */
#define GSL_EZERODIV  12  /* tried to divide by zero */
#define GSL_ETOL      13  /* user specified an invalid tolerance */
#define GSL_EUNDRFLW  14  /* underflow */
#define GSL_EOVRFLW   15  /* overflow (duh... ) */


/* just to make things slightly clearer */
enum {GSL_SUCCESS = 0, GSL_FAILURE = -1};

typedef void gsl_errhandler_t (const char *reason, const char *file, int line);

void gsl_error (const char *reason, const char *file, int line);
void (*gsl_set_error_handler (void (*new_handler) (const char *reason, const char *file, int line))) (const char *reason, const char *file, int line);

void gsl_empty_error_handler (const char *reason, const char * file, int line);
void gsl_message(const char * message, const char * file, int line);

#define GSL_ERROR(reason, gsl_errno) \
       do { \
       gsl_error (reason, __FILE__, __LINE__) ; \
       return gsl_errno ; \
       } while (0)

#define GSL_ERROR_RETURN(reason, gsl_errno, value) \
       do { \
       gsl_error (reason, __FILE__, __LINE__) ; \
       return value ; \
       } while (0)

#ifdef GSL_MESSAGING /* output all messages */
#define GSL_MESSAGE(message) \
       do { \
       gsl_message (message, __FILE__, __LINE__) ; \
       } while (0)
#else /* throw away message */
#define GSL_MESSAGE(message) do { } while(0)
#endif

#endif /* _GSL_ERRNO_H */
