noinst_LIBRARIES = libgslvector.a 

check_PROGRAMS = test test_static

include_HEADERS = 		\
	gsl_vector.h		\
	gsl_vector_char.h	\
	gsl_vector_complex.h	\
	gsl_vector_complex_double.h	\
	gsl_vector_complex_float.h	 \
	gsl_vector_complex_long_double.h \
	gsl_vector_double.h	\
	gsl_vector_float.h	\
	gsl_vector_int.h	\
	gsl_vector_long.h	\
	gsl_vector_long_double.h\
	gsl_vector_short.h	\
	gsl_vector_uchar.h	\
	gsl_vector_uint.h	\
	gsl_vector_ulong.h	\
	gsl_vector_ushort.h

INCLUDES = -I$(srcdir)/.. -I$(srcdir)/../err

TESTS = test test_static

test_LDADD = libgslvector.a ../err/libgslerr.a ../err/libgsltest.a
test_static_LDADD = libgslvector.a ../err/libgslerr.a ../err/libgsltest.a

test_DEPENDENCIES = libgslvector.a
test_static_DEPENDENCIES = libgslvector.a

test_SOURCES = test.c
test_static_SOURCES = test_static.c

CLEANFILES = test.txt test.dat

noinst_HEADERS = vector_source.c init_source.c file_source.c fprintf_block.c \
	test_complex_source.c test_source.c

libgslvector_a_SOURCES = init.c file.c vector.c
