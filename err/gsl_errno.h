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
#define GSL_EOVRFLW   15  /* overflow  */
#define GSL_ELOSS     16  /* loss of accuracy */


/* just to make things slightly clearer */
enum {GSL_SUCCESS = 0, GSL_FAILURE = -1};


typedef void gsl_errhandler_t (const char *reason, const char *file, int line);

void gsl_error (const char *reason, const char *file, int line);
void (*gsl_set_error_handler (void (*new_handler) (const char *reason, const char *file, int line))) (const char *reason, const char *file, int line);

void gsl_empty_error_handler (const char *reason, const char * file, int line);

void gsl_message(const char * message, const char * file, int line, unsigned int mask);


/* call error handler, returning error status in POSIX fashion */
#define GSL_ERROR(reason, gsl_errno) \
       do { \
       gsl_error (reason, __FILE__, __LINE__) ; \
       return gsl_errno ; \
       } while (0)

/* call error handler, returning a specified value [non-POSIX behaviour] */
#define GSL_ERROR_RETURN(reason, gsl_errno, value) \
       do { \
       gsl_error (reason, __FILE__, __LINE__) ; \
       return value ; \
       } while (0)



/* Provide a general messaging service for client use.
 * Messages can be selectively turned off at compile time
 * by defining an appropriate message mask. Client code
 * which uses the GSL_MESSAGE() macro must provide a mask
 * which is or'ed with the GSL_MESSAGE_MASK.
 *
 * The messaging service can be completely turned off
 * by defining GSL_MESSAGING_OFF.
 */
#ifndef GSL_MESSAGE_MASK
#define GSL_MESSAGE_MASK 0xffffffff  /* default all messages allowed */
#endif

/* Provide some symolic masks for client ease of use. */
enum {
  GSL_MESSAGE_MASK_A = 1,
  GSL_MESSAGE_MASK_B = 2,
  GSL_MESSAGE_MASK_C = 4,
  GSL_MESSAGE_MASK_D = 8,
  GSL_MESSAGE_MASK_E = 16,
  GSL_MESSAGE_MASK_F = 32,
  GSL_MESSAGE_MASK_G = 64,
  GSL_MESSAGE_MASK_H = 128
};

#ifdef GSL_MESSAGING_OFF        /* throw away messages */ 
#define GSL_MESSAGE(message, mask) do { } while(0)
#else                           /* output all messages */
#define GSL_MESSAGE(message, mask) \
       do { \
       gsl_message (message, __FILE__, __LINE__, mask) ; \
       } while (0)
#endif


/* GSL library code can occasionally generate warnings, which
 * are not intended to be fatal. Warnings can be turned off globally
 * by defining the preprocessor constant GSL_WARNINGS_OFF. This
 * turns off all warning level messages, but does not disable error
 * handling in any way or turn off error messages.
 *
 * GSL_WARNING() is not intended for use in client code.
 * Use GSL_MESSAGE() instead.
 */
#ifdef GSL_WARNINGS_OFF   /* throw away warnings */
#define GSL_WARNING(warning) do { } while(0)
#else                     /* output all warnings */
#define GSL_WARNING(warning) \
       do { \
       gsl_warning (warning, __FILE__, __LINE__) ; \
       } while (0)
#endif


#endif /* !_GSL_ERRNO_H */
