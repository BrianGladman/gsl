#include <config.h>
#include <gsl_matrix_complex.h>

static const gsl_complex zero = {0, 0} ;

#define BASE gsl_complex
#define SHORT complex
#define ZERO zero

#include "init_source.c"

