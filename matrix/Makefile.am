lib_LIBRARIES = libgslmatrix.a 

check_PROGRAMS = test test_trap

include_HEADERS = gsl_matrix.h gsl_matrix_float.h gsl_matrix_int.h

INCLUDES = -I$(srcdir)/.. -I$(srcdir)/../err -I$(srcdir)/../utils

TESTS = test test_trap

test_LDADD = libgslmatrix.a ../err/libgslerr.a ../err/libgsltest.a ../utils/libutils.a
test_trap_LDADD = libgslmatrix.a ../err/libgslerr.a ../err/libgsltest.a ../utils/libutils.a

test_DEPENDENCIES = libgslmatrix.a
test_trap_DEPENDENCIES = libgslmatrix.a

test_SOURCES = test.c
test_trap_SOURCES = test_trap.c

noinst_HEADERS = matrix_source.c init_source.c


libgslmatrix_a_SOURCES = \
source.h \
init.c \
matrix.c \
gsl_matrix.h \
init_int.c \
matrix_int.c \
gsl_matrix_int.h \
init_float.c \
matrix_float.c \
gsl_matrix_float.h



