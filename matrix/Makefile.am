noinst_LIBRARIES = libgslmatrix.a 

check_PROGRAMS = test test_static

include_HEADERS = 			\
	gsl_matrix.h 			\
	gsl_matrix_char.h		\
	gsl_matrix_complex_double.h	\
	gsl_matrix_complex_float.h	\
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

TESTS = test test_static 

test_LDADD = libgslmatrix.a ../vector/libgslvector.a ../err/libgslerr.a ../err/libgsltest.a
test_static_LDADD = libgslmatrix.a ../vector/libgslvector.a ../err/libgslerr.a ../err/libgsltest.a

test_DEPENDENCIES = libgslmatrix.a
test_SOURCES = test.c

test_static_DEPENDENCIES = libgslmatrix.a
test_static_SOURCES = test_static.c

CLEANFILES = test.txt test.dat

noinst_HEADERS = matrix_source.c init_source.c file_source.c rowcol_source.c \
	test_complex_source.c test_source.c

libgslmatrix_a_SOURCES = init.c matrix.c file.c rowcol.c


