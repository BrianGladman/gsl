pkglib_LIBRARIES = libgslsort.a

pkginclude_HEADERS = gsl_sort.h gsl_sort_vector.h gsl_sort_vector_char.h gsl_sort_vector_double.h gsl_sort_vector_float.h gsl_sort_vector_int.h gsl_sort_vector_long.h gsl_sort_vector_long_double.h gsl_sort_vector_short.h gsl_sort_vector_uchar.h gsl_sort_vector_uint.h gsl_sort_vector_ulong.h gsl_sort_vector_ushort.h

INCLUDES= -I$(top_builddir)

libgslsort_a_SOURCES = sort.c sortvec.c
noinst_HEADERS = sortvec_source.c sortind_source.c

TESTS = test

check_PROGRAMS = test

test_SOURCES = test.c
test_LDADD = libgslsort.a ../permutation/libgslpermutation.a ../vector/libgslvector.a ../block/libgslblock.a ../err/libgslerr.a ../err/libgsltest.a

