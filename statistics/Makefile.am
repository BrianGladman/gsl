## Process this file with automake to produce Makefile.in

lib_LIBRARIES = libgslstatistics.a
include_HEADERS = gsl_statistics.h

INCLUDES = -I$(srcdir)/.. -I$(srcdir)/../err

noinst_PROGRAMS = test

TESTS = test

test_SOURCES = test.c
test_LDADD = libgslstatistics.a ../err/libgsltest.a
test_DEPENDENCIES = libgslstatistics.a

libgslstatistics_a_SOURCES = \
mean.c \
variance.c variance.h \
pooled_variance.c \
minmax.c \
ttest.c \
imean.c \
ivariance.c ivariance.h \
ipooled_variance.c \
iminmax.c \
ittest.c




