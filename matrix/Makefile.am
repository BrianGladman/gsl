lib_LIBRARIES = libgslmatrix.a 

check_PROGRAMS = test test_int test_float test_complex test_trap

include_HEADERS = 			\
	gsl_matrix.h 			\
	gsl_matrix_char.h		\
	gsl_matrix_complex_double.h	\
	gsl_matrix_complex_long_double.h\
	gsl_matrix_double.h		\
	gsl_matrix_float.h		\
	gsl_matrix_int.h		\
	gsl_matrix_long.h		\
	gsl_matrix_long_double.h	\
	gsl_matrix_short.h		\
	gsl_matrix_uchar.h		\
	gsl_matrix_uint.h		\
	gsl_matrix_ulong.h		\
	gsl_matrix_ushort.h

INCLUDES = -I$(srcdir)/.. -I$(srcdir)/../vector -I$(srcdir)/../err

TESTS = test test_float test_int test_complex test_trap

test_LDADD = libgslmatrix.a ../vector/libgslvector.a ../err/libgslerr.a ../err/libgsltest.a
test_int_LDADD = libgslmatrix.a ../vector/libgslvector.a ../err/libgslerr.a ../err/libgsltest.a
test_float_LDADD = libgslmatrix.a ../vector/libgslvector.a ../err/libgslerr.a ../err/libgsltest.a
test_complex_LDADD = libgslmatrix.a ../vector/libgslvector.a ../err/libgslerr.a ../err/libgsltest.a
test_trap_LDADD = libgslmatrix.a ../vector/libgslvector.a ../err/libgslerr.a ../err/libgsltest.a

test_DEPENDENCIES = libgslmatrix.a
test_trap_DEPENDENCIES = libgslmatrix.a

test_SOURCES = test.c
test_int_SOURCES = test_int.c
test_float_SOURCES = test_float.c
test_complex_SOURCES = test_complex.c
test_trap_SOURCES = test_trap.c

CLEANFILES = test.txt test.dat

noinst_HEADERS = matrix_source.c init_source.c file_source.c source.h unsource.h

libgslmatrix_a_SOURCES = source.h init.c matrix.c file.c
