.PHONY: all
all:
	cd src; make all

.PHONY: gtest
gtest:
	cd gtest; make all

configure: configure.ac config.make.in
	autoconf
	autoheader
	-rm -fr autom4te.cache

.PHONY: htmldoc
htmldoc:
	rm -fr htmldoc
	doxygen Doxyfile

.PHONY: clean
clean:
	cd src; make clean
	cd gtest; make clean
	cd example; make clean

.PHONY: distclean
distclean: clean
	cd src; make distclean
	-rm -f config.status config.log config.hpp config.make
