## Process this file with automake to produce Makefile.in

lib_LIBRARIES = libgslstatistics.a
include_HEADERS = gsl_statistics.h gsl_statistics_int.h

INCLUDES = -I$(srcdir)/.. -I$(srcdir)/../err

noinst_PROGRAMS = test test_int
noinst_HEADERS = mean_source.c \
variance_source.c  \
absdev_source.c \
skew_source.c \
kurtosis_source.c \
pooled_variance_source.c \
minmax_source.c \
ttest_source.c \
sort_source.c \
median_source.c \
percentiles_source.c \
source.h

TESTS = test test_int

test_SOURCES = test.c test.h
test_LDADD = libgslstatistics.a ../err/libgsltest.a
test_DEPENDENCIES = libgslstatistics.a

test_int_SOURCES = test_int.c test.h
test_int_LDADD = libgslstatistics.a ../err/libgsltest.a
test_int_DEPENDENCIES = libgslstatistics.a

libgslstatistics_a_SOURCES = \
mean.c \
variance.c  \
absdev.c \
skew.c \
kurtosis.c \
pooled_variance.c \
minmax.c \
ttest.c \
sort.c \
median.c \
percentiles.c \
mean_int.c \
variance_int.c  \
absdev_int.c \
skew_int.c \
kurtosis_int.c \
pooled_variance_int.c \
minmax_int.c \
ttest_int.c \
sort_int.c \
median_int.c \
percentiles_int.c
