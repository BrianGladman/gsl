
header-links: remove-links
	HEADERLIST="$(top_srcdir)/gsl*.h $(top_srcdir)/*/gsl*.h"; \
	for h in $$HEADERLIST; do \
	  BASENAME=`basename $$h`; \
	  ln -s $$h $$BASENAME; \
	done

remove-links: 
	rm -rf $(full_top_builddir)/gsl/*.h


all: all-am header-links

clean: clean-am remove-links
