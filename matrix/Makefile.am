noinst_LTLIBRARIES = libgslmatrix.la 

check_PROGRAMS = test test_static

pkginclude_HEADERS = gsl_matrix.h gsl_matrix_char.h gsl_matrix_complex_double.h gsl_matrix_complex_float.h gsl_matrix_complex_long_double.h gsl_matrix_double.h gsl_matrix_float.h gsl_matrix_int.h gsl_matrix_long.h gsl_matrix_long_double.h gsl_matrix_short.h gsl_matrix_uchar.h gsl_matrix_uint.h gsl_matrix_ulong.h gsl_matrix_ushort.h

INCLUDES= -I$(top_builddir) -I$(top_srcdir)

TESTS = test test_static 

test_LDADD = libgslmatrix.la ../vector/libgslvector.la ../block/libgslblock.la ../err/libgslerr.la ../err/libgsltest.la ../utils/libutils.la
test_static_LDADD = libgslmatrix.la ../vector/libgslvector.la ../block/libgslblock.la ../err/libgslerr.la ../err/libgsltest.la ../utils/libutils.la

test_DEPENDENCIES = libgslmatrix.la
test_SOURCES = test.c

test_static_DEPENDENCIES = libgslmatrix.la
test_static_SOURCES = test_static.c

CLEANFILES = test.txt test.dat

noinst_HEADERS = matrix_source.c init_source.c file_source.c rowcol_source.c swap_source.c copy_source.c test_complex_source.c test_source.c minmax_source.c prop_source.c

libgslmatrix_la_SOURCES = init.c matrix.c file.c rowcol.c swap.c copy.c minmax.c prop.c


