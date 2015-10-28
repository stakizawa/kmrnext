.PHONY: all
all:
	cd src; make all

configure: configure.ac config.make.in
	autoconf
	autoheader
	-rm -fr autom4te.cache

.PHONY: clean
clean:
	cd src; make clean

.PHONY: distclean
distclean:
	cd src; make distclean
	-rm -f config.status config.log config.h config.make
