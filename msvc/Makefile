all:
	perl am2ms.pl ../Makefile.am ../*/Makefile.am
	#perl -i -p -e 's/\012/\015\012/' *.dsp *.dsw *.bat
	mv *.dsw ..
	mv *.bat ..
	perl gensetup.pl setup.iss.in > setup.iss
	perl mkexports.pl libgsl.nm.dat > libgsl.def
	perl mkexports.pl libgslcblas.nm.dat > libgslcblas.def
        
update:
	nm ../.libs/libgsl.a > libgsl.nm.dat
	wc -l libgsl.nm.dat
	nm ../cblas/.libs/libgslcblas.a > libgslcblas.nm.dat
	wc -l libgslcblas.nm.dat

clean:
	rm -f *.dsp *.dsw *.bat *.plg *.opt *.ncb

