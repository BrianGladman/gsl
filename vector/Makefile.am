lib_LIBRARIES = libgslvector.a 

check_PROGRAMS = test test_float test_int test_complex test_trap 

include_HEADERS = 		\
	gsl_vector.h		\
	gsl_vector_char.h	\
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

TESTS = test test_float test_int test_complex test_trap

test_LDADD = libgslvector.a ../err/libgslerr.a ../err/libgsltest.a
test_float_LDADD = libgslvector.a ../err/libgslerr.a ../err/libgsltest.a
test_int_LDADD = libgslvector.a ../err/libgslerr.a ../err/libgsltest.a
test_complex_LDADD = libgslvector.a ../err/libgslerr.a ../err/libgsltest.a
test_trap_LDADD = libgslvector.a ../err/libgslerr.a ../err/libgsltest.a

test_DEPENDENCIES = libgslvector.a
test_trap_DEPENDENCIES = libgslvector.a

test_SOURCES = test.c
test_float_SOURCES = test_float.c
test_int_SOURCES = test_int.c
test_int_SOURCES = test_complex.c
test_trap_SOURCES = test_trap.c

CLEANFILES = test.txt test.dat

noinst_HEADERS = vector_source.c init_source.c file_source.c fprintf_block.c

libgslvector_a_SOURCES = init.c file.c vector.c
