all:
	perl am2ms.pl ../Makefile.am ../*/Makefile.am
	perl -i -p -e 's/\012/\015\012/' *.dsp *.dsw *.bat
	mv *.dsw ..
	mv *.bat ..

clean:
	rm -f *.dsp *.dsw *.bat *.plg *.opt *.ncb

zip:
	rm -f GSL_MSVC_BUILD.zip
	(cd .. ; zip msvc/GSL_MSVC_BUILD.zip msvc/*.{txt,bat,dsw,dsp})
