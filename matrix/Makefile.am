lib_LIBRARIES = libgslmatrix.a 

check_PROGRAMS = test test_int test_float test_trap

include_HEADERS = gsl_matrix.h gsl_matrix_float.h gsl_matrix_int.h

INCLUDES = -I$(srcdir)/.. -I$(srcdir)/../vector -I$(srcdir)/../err -I$(srcdir)/../utils

TESTS = test test_float test_int test_trap

test_LDADD = libgslmatrix.a ../vector/libgslvector.a ../err/libgslerr.a ../err/libgsltest.a ../utils/libutils.a
test_int_LDADD = libgslmatrix.a ../vector/libgslvector.a ../err/libgslerr.a ../err/libgsltest.a ../utils/libutils.a
test_float_LDADD = libgslmatrix.a ../vector/libgslvector.a ../err/libgslerr.a ../err/libgsltest.a ../utils/libutils.a
test_trap_LDADD = libgslmatrix.a ../vector/libgslvector.a ../err/libgslerr.a ../err/libgsltest.a ../utils/libutils.a

test_DEPENDENCIES = libgslmatrix.a
test_trap_DEPENDENCIES = libgslmatrix.a

test_SOURCES = test.c
test_int_SOURCES = test_int.c
test_float_SOURCES = test_float.c
test_trap_SOURCES = test_trap.c

noinst_HEADERS = matrix_source.c init_source.c file_source.c


libgslmatrix_a_SOURCES = \
source.h \
init.c \
matrix.c \
file.c \
gsl_matrix.h \
init_int.c \
matrix_int.c \
file_int.c \
gsl_matrix_int.h \
init_float.c \
matrix_float.c \
file_float.c \
gsl_matrix_float.h



