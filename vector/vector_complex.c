#include <config.h>
#include <gsl_vector_complex.h>

static gsl_complex zero = {0, 0} ;

#define BASE gsl_complex
#define SHORT complex
#define ZERO zero
#include "vector_source.c"

