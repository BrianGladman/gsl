# Generated automatically from Makefile.in by configure.
# Makefile.in generated automatically by automake 1.2 from Makefile.am

# Copyright (C) 1994, 1995, 1996, 1997 Free Software Foundation, Inc.
# This Makefile.in is free software; the Free Software Foundation
# gives unlimited permission to copy, distribute and modify it.


SHELL = /bin/sh

srcdir = .
top_srcdir = ..
prefix = /usr/local
exec_prefix = ${prefix}

bindir = ${exec_prefix}/bin
sbindir = ${exec_prefix}/sbin
libexecdir = ${exec_prefix}/libexec
datadir = ${prefix}/share
sysconfdir = ${prefix}/etc
sharedstatedir = ${prefix}/com
localstatedir = ${prefix}/var
libdir = ${exec_prefix}/lib
infodir = ${prefix}/info
mandir = ${prefix}/man
includedir = ${prefix}/include
oldincludedir = /usr/include

pkgdatadir = $(datadir)/gsl
pkglibdir = $(libdir)/gsl
pkgincludedir = $(includedir)/gsl

top_builddir = ..

ACLOCAL = aclocal
AUTOCONF = autoconf
AUTOMAKE = automake
AUTOHEADER = autoheader

INSTALL = /usr/bin/install -c
INSTALL_PROGRAM = ${INSTALL}
INSTALL_DATA = ${INSTALL} -m 644
INSTALL_SCRIPT = ${INSTALL_PROGRAM}
transform = s,x,x,

NORMAL_INSTALL = true
PRE_INSTALL = true
POST_INSTALL = true
NORMAL_UNINSTALL = true
PRE_UNINSTALL = true
POST_UNINSTALL = true
host_alias = i586-pc-linux-gnu
host_triplet = i586-pc-linux-gnu
CC = gcc
CPP = gcc -E
LN_S = ln -s
MAKEINFO = makeinfo
PACKAGE = gsl
RANLIB = ranlib
VERSION = 0.3a-interim

lib_LIBRARIES = libgslspecfunc.a 

# noinst_PROGRAMS = 

include_HEADERS = \
	gsl_sf.h 	    gsl_sf_chebyshev.h  gsl_sf_gamma.h      gsl_sf_pow_int.h \
	gsl_sf_airy.h       gsl_sf_dilog.h      gsl_sf_laguerre.h   gsl_sf_trig.h \
	gsl_sf_bessel.h     gsl_sf_erf.h        gsl_sf_log.h \
	gsl_sf_legendre.h

noinst_HEADERS = bessel_amp_phase.h

INCLUDES = -I$(srcdir)/../err -I$(srcdir)/../integration

TESTS = test_trap test_radix2 test

#test_LDADD = libgslfft.a ../err/libgslerr.a ../err/libgsltest.a ../utils/libutils.a
#test_trap_LDADD = libgslfft.a ../err/libgslerr.a ../err/libgsltest.a ../utils/libutils.a
#test_radix2_LDADD = libgslfft.a ../err/libgslerr.a ../err/libgsltest.a ../utils/libutils.a
#errs_LDADD = libgslfft.a ../err/libgslerr.a ../err/libgsltest.a ../utils/libutils.a
#benchmark_LDADD = libgslfft.a ../err/libgslerr.a ../err/libgsltest.a ../utils/libutils.a

#test_DEPENDENCIES = libgslfft.a
#test_trap_DEPENDENCIES = libgslfft.a
#test_radix2_DEPENDENCIES = libgslfft.a
#benchmark_DEPENDENCIES = libgslfft.a

#test_SOURCES = test.c
#test_trap_SOURCES = test_trap.c
#test_radix2_SOURCES = test_radix2.c
#benchmark_SOURCES = benchmark.c 

libgslspecfunc_a_SOURCES = 	\
	airy.c			\
	bessel.c		\
	bessel_I0.c		\
	bessel_I1.c		\
	bessel_J0.c		\
	bessel_J1.c		\
	bessel_K0.c		\
	bessel_K1.c		\
	bessel_Y0.c		\
	bessel_Y1.c		\
	bessel_amp_phase.c	\
	chebyshev.c		\
	dilog.c			\
	erfc.c			\
	gamma.c			\
	laguerre.c		\
	legendre.c		\
	log.c			\
	pow_int.c		\
	sig2lamp.c		\
	t-erf.c			\
	trig.c
mkinstalldirs = $(SHELL) $(top_srcdir)/mkinstalldirs
CONFIG_HEADER = ../config.h
CONFIG_CLEAN_FILES = 
LIBRARIES =  $(lib_LIBRARIES)


DEFS = -DHAVE_CONFIG_H -I. -I$(srcdir) -I..
CPPFLAGS = 
LDFLAGS = 
LIBS = -lm 
libgslspecfunc_a_LIBADD = 
libgslspecfunc_a_OBJECTS =  airy.o bessel.o bessel_I0.o bessel_I1.o \
bessel_J0.o bessel_J1.o bessel_K0.o bessel_K1.o bessel_Y0.o bessel_Y1.o \
bessel_amp_phase.o chebyshev.o dilog.o erfc.o gamma.o laguerre.o \
legendre.o log.o pow_int.o sig2lamp.o t-erf.o trig.o
AR = ar
CFLAGS = -g -O2
COMPILE = $(CC) $(DEFS) $(INCLUDES) $(CPPFLAGS) $(CFLAGS)
LINK = $(CC) $(CFLAGS) $(LDFLAGS) -o $@
HEADERS =  $(include_HEADERS) $(noinst_HEADERS)

DIST_COMMON =  Makefile.am Makefile.in


DISTFILES = $(DIST_COMMON) $(SOURCES) $(HEADERS) $(TEXINFOS) $(EXTRA_DIST)

TAR = tar
GZIP = --best
DEP_FILES =  .deps/airy.P .deps/bessel.P .deps/bessel_I0.P \
.deps/bessel_I1.P .deps/bessel_J0.P .deps/bessel_J1.P .deps/bessel_K0.P \
.deps/bessel_K1.P .deps/bessel_Y0.P .deps/bessel_Y1.P \
.deps/bessel_amp_phase.P .deps/chebyshev.P .deps/dilog.P .deps/erfc.P \
.deps/gamma.P .deps/laguerre.P .deps/legendre.P .deps/log.P \
.deps/pow_int.P .deps/sig2lamp.P .deps/t-erf.P .deps/trig.P
SOURCES = $(libgslspecfunc_a_SOURCES)
OBJECTS = $(libgslspecfunc_a_OBJECTS)

default: all

.SUFFIXES:
.SUFFIXES: .c .o
$(srcdir)/Makefile.in: Makefile.am $(top_srcdir)/configure.in $(ACLOCAL_M4) 
	cd $(top_srcdir) && $(AUTOMAKE) --gnu specfunc/Makefile

Makefile: $(srcdir)/Makefile.in $(top_builddir)/config.status $(BUILT_SOURCES)
	cd $(top_builddir) \
	  && CONFIG_FILES=$(subdir)/$@ CONFIG_HEADERS= $(SHELL) ./config.status


mostlyclean-libLIBRARIES:

clean-libLIBRARIES:
	test -z "$(lib_LIBRARIES)" || rm -f $(lib_LIBRARIES)

distclean-libLIBRARIES:

maintainer-clean-libLIBRARIES:

install-libLIBRARIES: $(lib_LIBRARIES)
	@$(NORMAL_INSTALL)
	$(mkinstalldirs) $(libdir)
	list='$(lib_LIBRARIES)'; for p in $$list; do \
	  if test -f $$p; then \
	    echo " $(INSTALL_DATA) $$p $(libdir)/$$p"; \
	    $(INSTALL_DATA) $$p $(libdir)/$$p; \
	  else :; fi; \
	done
	@$(POST_INSTALL)
	@list='$(lib_LIBRARIES)'; for p in $$list; do \
	  if test -f $$p; then \
	    echo " $(RANLIB) $(libdir)/$$p"; \
	    $(RANLIB) $(libdir)/$$p; \
	  else :; fi; \
	done

uninstall-libLIBRARIES:
	$(NORMAL_UNINSTALL)
	list='$(lib_LIBRARIES)'; for p in $$list; do \
	  rm -f $(libdir)/$$p; \
	done

.c.o:
	$(COMPILE) -c $<

mostlyclean-compile:
	rm -f *.o core

clean-compile:

distclean-compile:
	rm -f *.tab.c

maintainer-clean-compile:

libgslspecfunc.a: $(libgslspecfunc_a_OBJECTS) $(libgslspecfunc_a_DEPENDENCIES)
	rm -f libgslspecfunc.a
	$(AR) cru libgslspecfunc.a $(libgslspecfunc_a_OBJECTS) $(libgslspecfunc_a_LIBADD)
	$(RANLIB) libgslspecfunc.a

install-includeHEADERS: $(include_HEADERS)
	@$(NORMAL_INSTALL)
	$(mkinstalldirs) $(includedir)
	@list='$(include_HEADERS)'; for p in $$list; do \
	  echo " $(INSTALL_DATA) $(srcdir)/$$p $(includedir)/$$p"; \
	  $(INSTALL_DATA) $(srcdir)/$$p $(includedir)/$$p; \
	done

uninstall-includeHEADERS:
	$(NORMAL_UNINSTALL)
	list='$(include_HEADERS)'; for p in $$list; do \
	  rm -f $(includedir)/$$p; \
	done

tags: TAGS

ID: $(HEADERS) $(SOURCES)
	here=`pwd` && cd $(srcdir) && mkid -f$$here/ID $(SOURCES) $(HEADERS)

TAGS:  $(HEADERS) $(SOURCES)  $(TAGS_DEPENDENCIES)
	tags=; \
	here=`pwd`; \
	test -z "$(ETAGS_ARGS)$(SOURCES)$(HEADERS)$$tags" \
	  || (cd $(srcdir) && etags $(ETAGS_ARGS) $$tags  $(SOURCES) $(HEADERS) -o $$here/TAGS)

mostlyclean-tags:

clean-tags:

distclean-tags:
	rm -f TAGS ID

maintainer-clean-tags:

distdir = $(top_builddir)/$(PACKAGE)-$(VERSION)/$(subdir)

subdir = specfunc

distdir: $(DISTFILES)
	here=`cd $(top_builddir) && pwd`; \
	top_distdir=`cd $(top_distdir) && pwd`; \
	cd $(top_srcdir) \
	  && $(AUTOMAKE) --include-deps --build-dir=$$here --srcdir-name=$(top_srcdir) --output-dir=$$top_distdir --gnu specfunc/Makefile
	@for file in $(DISTFILES); do \
	  d=$(srcdir); \
	  test -f $(distdir)/$$file \
	  || ln $$d/$$file $(distdir)/$$file 2> /dev/null \
	  || cp -p $$d/$$file $(distdir)/$$file; \
	done

MKDEP = gcc -M $(DEFS) $(INCLUDES) $(CPPFLAGS) $(CFLAGS)

DEPS_MAGIC := $(shell mkdir .deps > /dev/null 2>&1 || :)
-include .deps/.P
.deps/.P: $(BUILT_SOURCES)
	echo > $@

-include $(DEP_FILES)

mostlyclean-depend:

clean-depend:

distclean-depend:

maintainer-clean-depend:
	rm -rf .deps

.deps/%.P: %.c
	@echo "Computing dependencies for $<..."
	@o='o'; \
	test -n "$o" && o='$$o'; \
	$(MKDEP) $< >$@.tmp \
	  && sed "s,^\(.*\)\.o:,\1.$$o \1.l$$o $@:," < $@.tmp > $@ \
	  && rm -f $@.tmp
check-TESTS: $(TESTS)
	@failed=0; all=0; \
	srcdir=$(srcdir); export srcdir; \
	for tst in $(TESTS); do \
	  if test -f $$tst; then dir=.; \
	  else dir="$(srcdir)"; fi; \
	  if $(TESTS_ENVIRONMENT) $$dir/$$tst; then \
	    all=`expr $$all + 1`; \
	    echo "PASS: $$tst"; \
	  elif test $$? -ne 77; then \
	    all=`expr $$all + 1`; \
	    failed=`expr $$failed + 1`; \
	    echo "FAIL: $$tst"; \
	  fi; \
	done; \
	if test "$$failed" -eq 0; then \
	  banner="All $$all tests passed"; \
	else \
	  banner="$$failed of $$all tests failed"; \
	fi; \
	dashes=`echo "$$banner" | sed s/./=/g`; \
	echo "$$dashes"; \
	echo "$$banner"; \
	echo "$$dashes"; \
	test "$$failed" -eq 0
info:
dvi:
check: all
	$(MAKE) check-TESTS
installcheck:
install-exec: install-libLIBRARIES
	@$(NORMAL_INSTALL)

install-data: install-includeHEADERS
	@$(NORMAL_INSTALL)

install: install-exec install-data all
	@:

uninstall: uninstall-libLIBRARIES uninstall-includeHEADERS

all: Makefile $(LIBRARIES) $(HEADERS)

install-strip:
	$(MAKE) INSTALL_PROGRAM='$(INSTALL_PROGRAM) -s' INSTALL_SCRIPT='$(INSTALL_PROGRAM)' install
installdirs:
	$(mkinstalldirs)  $(libdir) $(includedir)


mostlyclean-generic:
	test -z "$(MOSTLYCLEANFILES)" || rm -f $(MOSTLYCLEANFILES)

clean-generic:
	test -z "$(CLEANFILES)" || rm -f $(CLEANFILES)

distclean-generic:
	rm -f Makefile $(DISTCLEANFILES)
	rm -f config.cache config.log stamp-h stamp-h[0-9]*
	test -z "$(CONFIG_CLEAN_FILES)" || rm -f $(CONFIG_CLEAN_FILES)

maintainer-clean-generic:
	test -z "$(MAINTAINERCLEANFILES)" || rm -f $(MAINTAINERCLEANFILES)
	test -z "$(BUILT_SOURCES)" || rm -f $(BUILT_SOURCES)
mostlyclean:  mostlyclean-libLIBRARIES mostlyclean-compile \
		mostlyclean-tags mostlyclean-depend mostlyclean-generic

clean:  clean-libLIBRARIES clean-compile clean-tags clean-depend \
		clean-generic mostlyclean

distclean:  distclean-libLIBRARIES distclean-compile distclean-tags \
		distclean-depend distclean-generic clean
	rm -f config.status

maintainer-clean:  maintainer-clean-libLIBRARIES \
		maintainer-clean-compile maintainer-clean-tags \
		maintainer-clean-depend maintainer-clean-generic \
		distclean
	@echo "This command is intended for maintainers to use;"
	@echo "it deletes files that may require special tools to rebuild."

.PHONY: default mostlyclean-libLIBRARIES distclean-libLIBRARIES \
clean-libLIBRARIES maintainer-clean-libLIBRARIES uninstall-libLIBRARIES \
install-libLIBRARIES mostlyclean-compile distclean-compile \
clean-compile maintainer-clean-compile uninstall-includeHEADERS \
install-includeHEADERS tags mostlyclean-tags distclean-tags clean-tags \
maintainer-clean-tags distdir mostlyclean-depend distclean-depend \
clean-depend maintainer-clean-depend check-TESTS info dvi installcheck \
install-exec install-data install uninstall all installdirs \
mostlyclean-generic distclean-generic clean-generic \
maintainer-clean-generic clean mostlyclean distclean maintainer-clean


# Tell versions [3.59,3.63) of GNU make to not export all variables.
# Otherwise a system limit (for SysV at least) may be exceeded.
.NOEXPORT:
