all:
	for n in GSLLIB GSLLIBMT GSLLIBMD GSLDLL ; do  \
	( cd .. ; mkdir $$n ; cd $$n ; \
	  perl ../msvc/am2ms.pl $$n ../Makefile.am ../*/Makefile.am ; \
	  mv GSL.dsw ../$$n.dsw ; \
	  mv GSLTESTS.dsw ../$${n}TESTS.dsw ; \
	 ) \
	done
	#mv *.dsw ..
	#mv *.bat ..
	perl gensetup.pl setup.iss.in > setup.iss
	perl mkexports.pl libgsl.nm.dat > ../GSLDLL/libgsl.def
	perl mkexports.pl libgslcblas.nm.dat > ../GSLDLL/libgslcblas.def
	rm libgsl.nm.dat libgslcblas.nm.dat
        
update:
	nm ../.libs/libgsl.a > libgsl.nm.dat
	wc -l libgsl.nm.dat
	nm ../cblas/.libs/libgslcblas.a > libgslcblas.nm.dat
	wc -l libgslcblas.nm.dat

clean:
	rm -f *.dsp *.dsw *.bat *.plg *.opt *.ncb
	rm -f libgsl.nm.dat libgslcblas.nm.dat

