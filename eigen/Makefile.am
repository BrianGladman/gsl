noinst_LTLIBRARIES = libgsleigen.la 

check_PROGRAMS = test

pkginclude_HEADERS = \
	gsl_eigen.h

# noinst_HEADERS = tridiag.h

INCLUDES= -I$(top_builddir)

TESTS = test

test_LDADD = libgsleigen.la 		\
	../sys/libgslsys.la		\
	../err/libgsltest.la		\
	../linalg/libgsllinalg.la        \
	../matrix/libgslmatrix.la	\
	../vector/libgslvector.la	\
	../block/libgslblock.la	        \
 	../err/libgslerr.la              \
        ../utils/libutils.la

test_DEPENDENCIES = libgsleigen.la
test_SOURCES = test.c

libgsleigen_la_SOURCES = eigen_sort.c jacobi.c
