## Process this file with automake to produce Makefile.in

pkglib_LIBRARIES = libgslstatistics.a
pkginclude_HEADERS = gsl_statistics.h \
 gsl_statistics_char.h gsl_statistics_double.h gsl_statistics_float.h \
 gsl_statistics_int.h gsl_statistics_long.h gsl_statistics_long_double.h \
 gsl_statistics_short.h gsl_statistics_uchar.h gsl_statistics_uint.h \
 gsl_statistics_ulong.h gsl_statistics_ushort.h

INCLUDES = -I$(srcdir)/.. -I$(srcdir)/../err

noinst_HEADERS = mean_source.c variance_source.c absdev_source.c \
 skew_source.c kurtosis_source.c lag1_source.c p_variance_source.c \
 minmax_source.c ttest_source.c sort_source.c median_source.c \
 quantiles_source.c source.h \
 test_float_source.c test_int_source.c

check_PROGRAMS = test test_nist
TESTS = test test_nist

test_SOURCES = test.c
test_LDADD = libgslstatistics.a ../err/libgsltest.a
test_DEPENDENCIES = libgslstatistics.a

test_nist_SOURCES = test_nist.c
test_nist_LDADD = libgslstatistics.a ../err/libgsltest.a
test_nist_DEPENDENCIES = libgslstatistics.a

libgslstatistics_a_SOURCES =  mean.c variance.c absdev.c skew.c \
 kurtosis.c lag1.c p_variance.c minmax.c ttest.c sort.c median.c \
 quantiles.c
