pkglib_LIBRARIES = libgsleigen.a 

check_PROGRAMS = test

pkginclude_HEADERS = \
	gsl_eigen.h

# noinst_HEADERS = tridiag.h

INCLUDES = \
	-I$(srcdir)/..		\
	-I$(srcdir)/../err	\
	-I$(srcdir)/../linalg	\
	-I$(srcdir)/../permutation \
	-I$(srcdir)/../matrix	\
	-I$(srcdir)/../vector   \
	-I$(srcdir)/../block

TESTS = test

test_LDADD = libgsleigen.a 		\
	../sys/libgslsys.a		\
	../err/libgsltest.a		\
	../linalg/libgsllinalg.a        \
	../matrix/libgslmatrix.a	\
	../vector/libgslvector.a	\
	../block/libgslblock.a	        \
 	../err/libgslerr.a              \
        ../utils/libutils.a

test_DEPENDENCIES = libgsleigen.a
test_SOURCES = test.c

libgsleigen_a_SOURCES = eigen_sort.c jacobi.c
