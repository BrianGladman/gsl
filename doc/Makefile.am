## Process this file with automake to produce Makefile.in

SUBDIRS = . examples

#info_TEXINFOS = gsl-ref.texi
noinst_TEXINFOS = gsl-design.texi

man_MANS = gsl.3 gsl-config.1 gsl-randist.1 gsl-histogram.1

figures =                              \
  images/bspline.png                   \
  images/cheb.png                      \
  images/dwt.png                       \
  images/fft-complex-radix2.png        \
  images/fft-real-mixedradix.png       \
  images/fit-exp.png                   \
  images/fit-wlinear2.png              \
  images/fit-wlinear.png               \
  images/histogram.png                 \
  images/interp2d.png                  \
  images/interp_compare.png            \
  images/interp.png                    \
  images/interpp.png                   \
  images/min-interval.png              \
  images/multilarge.png                \
  images/multimin.png                  \
  images/nlfit2b.png                   \
  images/nlfit2.png                    \
  images/nlfit3.png                    \
  images/ntuple.png                    \
  images/ode-vdp.png                   \
  images/qrng.png                      \
  images/rand-bernoulli.png            \
  images/rand-beta.png                 \
  images/rand-binomial.png             \
  images/rand-bivariate-gaussian.png   \
  images/rand-cauchy.png               \
  images/rand-chisq.png                \
  images/rand-erlang.png               \
  images/rand-exponential.png          \
  images/rand-exppow.png               \
  images/rand-fdist.png                \
  images/rand-flat.png                 \
  images/rand-gamma.png                \
  images/rand-gaussian.png             \
  images/rand-gaussian-tail.png        \
  images/rand-geometric.png            \
  images/rand-gumbel1.png              \
  images/rand-gumbel2.png              \
  images/rand-hypergeometric.png       \
  images/rand-landau.png               \
  images/rand-laplace.png              \
  images/rand-levy.png                 \
  images/rand-levyskew.png             \
  images/rand-logarithmic.png          \
  images/rand-logistic.png             \
  images/rand-lognormal.png            \
  images/rand-nbinomial.png            \
  images/random-walk.png               \
  images/rand-pareto.png               \
  images/rand-pascal.png               \
  images/rand-poisson.png              \
  images/rand-rayleigh.png             \
  images/rand-rayleigh-tail.png        \
  images/rand-tdist.png                \
  images/rand-weibull.png              \
  images/regularized2.png              \
  images/regularized.png               \
  images/robust.png                    \
  images/siman-12-cities.png           \
  images/siman-energy.png              \
  images/siman-final-route.png         \
  images/siman-initial-route.png       \
  images/siman-test.png                \
  images/sparse_poisson.png

noinst_DATA = $(figures)

EXTRA_DIST = $(man_MANS) $(noinst_DATA) gsl-design.texi gsl-ref.info

SPHINX_OPTS        =
SPHINX_BUILD       = sphinx-build
SPHINX_SOURCEDIR   = .
SPHINX_BUILDDIR    = _build

# targets supported by sphinx
SPHINX_TARGETS = html dirhtml singlehtml pickle json htmlhelp qthelp devhelp epub latex latexpdf latexpdfja text man texinfo gettext changes xml pseudoxml linkcheck doctest coverage

$(SPHINX_TARGETS) :
	@$(SPHINX_BUILD) -M $@ "$(SPHINX_SOURCEDIR)" "$(SPHINX_BUILDDIR)" $(SPHINX_OPTS) $(O)

# special case info target
info:
	@$(SPHINX_BUILD) -M $@ "$(SPHINX_SOURCEDIR)" "$(SPHINX_BUILDDIR)" $(SPHINX_OPTS) $(O)
	cp "$(SPHINX_BUILDDIR)/texinfo/gsl-ref.info" .

install-data-local: gsl-ref.info
	${INSTALL} gsl-ref.info ${infodir}

clean-local:
	rm -rf "$(SPHINX_BUILDDIR)"
