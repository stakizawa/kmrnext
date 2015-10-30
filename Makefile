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

.PHONY: clean
clean:
	cd src; make clean
	cd gtest; make clean

.PHONY: distclean
distclean: clean
	cd src; make distclean
	-rm -f config.status config.log config.h config.make
