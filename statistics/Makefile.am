## Process this file with automake to produce Makefile.in

lib_LIBRARIES = libgslstatistics.a
include_HEADERS = gsl_statistics.h gsl_statistics_int.h

INCLUDES = -I$(srcdir)/.. -I$(srcdir)/../err

noinst_PROGRAMS = test test_int

TESTS = test test_int

test_SOURCES = test.c
test_LDADD = libgslstatistics.a ../err/libgsltest.a
test_DEPENDENCIES = libgslstatistics.a

test_int_SOURCES = test_int.c
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
int_mean.c \
int_variance.c  \
int_absdev.c \
int_skew.c \
int_kurtosis.c \
int_pooled_variance.c \
int_minmax.c \
int_ttest.c \
int_sort.c \
int_median.c \
int_percentiles.c

int:
	perl convert_double_to_int.pl gsl_statistics.h > gsl_statistics_int.h 
	perl convert_double_to_int.pl mean.c > int_mean.c
	perl convert_double_to_int.pl variance.c > int_variance.c
	perl convert_double_to_int.pl absdev.c > int_absdev.c
	perl convert_double_to_int.pl minmax.c > int_minmax.c
	perl convert_double_to_int.pl pooled_variance.c > int_pooled_variance.c
	perl convert_double_to_int.pl ttest.c > int_ttest.c
	perl convert_double_to_int.pl skew.c > int_skew.c
	perl convert_double_to_int.pl kurtosis.c > int_kurtosis.c
	perl convert_double_to_int.pl sort.c > int_sort.c
	perl convert_double_to_int.pl median.c > int_median.c
	perl convert_double_to_int.pl percentiles.c > int_percentiles.c

