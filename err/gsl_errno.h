#include <errno.h>

#define GSL_EDOM   1		/* domain error for input, e.g. asin(10),sqrt(-1) */
#define GSL_ERANGE 2		/* range error for output, e.g. exp(1e100) */
#define GSL_EFAULT 3
#define GSL_EINVAL 4
#define GSL_EFAILED 5		/* generic failure */
#define GSL_EFACTOR 6		/* factorization failed */
#define GSL_ESANITY 7		/* sanity check failed -- this shouldn't happen */
#define GSL_ENOMEM 8		/* malloc failed */

extern int gsl_errno;

typedef void gsl_errhandler_t (const char *reason, const char *file, int line);

void gsl_error (const char *reason, const char *file, int line);
void (*gsl_set_error_handler (void (*new_handler) (const char *reason, const char *file, int line))) (const char *reason, const char *file, int line);

void gsl_no_error_handler (const char *reason, const char * file, int line) ;

#define GSL_ERROR(reason, errno) \
       do { \
       gsl_errno = errno ; \
       gsl_error (reason, __FILE__, __LINE__) ; \
       return -1 ; \
       } while (0) ;

#define GSL_ERROR_RETURN(reason, errno, value) \
       do { \
       gsl_errno = errno ; \
       gsl_error (reason, __FILE__, __LINE__) ; \
       return value ; \
       } while (0) ;

