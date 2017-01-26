include ../config.make

HEADERS = kmrnext.hpp ckmrnext.h kmrnextf.f90 \
	  data_element.hpp data_store.hpp \
	  data_store_serial.hpp data_store_kmr.hpp \
	  util.hpp
LIB     = libkmrnext.a
OBJS    = kmrnext.o basic_data.o data_store.o data_element.o ckmrnext.o \
	  api_fortran.o kmrnextf.o util.o
TESTS   = test00-basic test01-mng-ds test02-clang test03-fortran test04-file \
	  test05-multi

ifeq ($(TARGET), KMR)
    DEPLIB += $(KMRDIR)/lib/libkmr.a
endif

.Phony: all
all: $(LIB)

.Phony: tests
tests: $(TESTS)

libkmrnext.a: .update_header $(HEADERS) $(OBJS)
	ar csr $@ $(OBJS)

.SUFFIXES: .f90
%.o: %.mod

.cpp.o:
ifeq ($(TARGET), KMR)
	$(CXX) $(DEFS) $(CXXFLAGS) -I$(KMRDIR)/include -c $<
else
	$(CXX) $(DEFS) $(CXXFLAGS) -c $<
endif

.c.o:
ifeq ($(TARGET), KMR)
	$(CC) $(DEFS) $(CFLAGS) -I$(KMRDIR)/include -c $<
else
	$(CC) $(DEFS) $(CFLAGS) -c $<
endif

.f90.o:
ifeq ($(TARGET), KMR)
	$(FC) $(DEFS) -DBACKEND_KMR $(FCFLAGS) -I$(KMRDIR)/include -c $<
else
	$(FC) $(DEFS) $(FCFLAGS) -c $<
endif

test00-basic: test00-basic.o $(LIB)
	$(CXX) $(LDFLAGS) -o $@ $^ $(DEPLIB)

test01-mng-ds: test01-mng-ds.o $(LIB)
	$(CXX) $(LDFLAGS) -o $@ $^ $(DEPLIB)

test02-clang: test02-clang.o $(LIB)
	$(CC) $(LDFLAGS) -o $@ $^ $(DEPLIB) $(CXXLIB)

test03-fortran: test03-fortran.o test03-fortran-helper.o $(LIB)
	$(FC) $(LDFLAGS) -o $@ $^ $(DEPLIB) $(CXXLIB)

test04-file: test04-file.o $(LIB)
	$(CXX) $(LDFLAGS) -o $@ $^ $(DEPLIB)

test05-multi: test05-multi.o $(LIB)
	$(CXX) $(LDFLAGS) -o $@ $^ $(DEPLIB)

.update_header:
	cp kmrnext.hpp _kmrnext.hpp
	sed -e "s/#define BACKEND_[A-Z]* 1/#define BACKEND_$(TARGET) 1/" \
	< _kmrnext.hpp > kmrnext.hpp
	cp ckmrnext.h _ckmrnext.h
	sed -e "s/#define BACKEND_[A-Z]* 1/#define BACKEND_$(TARGET) 1/" \
	< _ckmrnext.h > ckmrnext.h
	cp kmrnextf.f90 _kmrnextf.f90
	sed -e "s/#define BACKEND_[A-Z]* 1/#define BACKEND_$(TARGET) 1/" \
	< _kmrnextf.f90 > kmrnextf.f90
	touch .update_header
	rm -f _kmrnext.hpp _ckmrnext.h _kmrnextf.f90

.PHONY: install
install: all
	install -d $(LIBDIR)
	install -d $(INCDIR)
	install -m 644 kmrnext.hpp $(INCDIR)
	install -m 644 ckmrnext.h $(INCDIR)
	install -m 644 kmrnextf.mod $(INCDIR)
	install -m 644 $(LIB) $(LIBDIR)

.PHONY: clean
clean:
	-rm -rf $(OBJS) $(TESTS) *.o *.d *.mod

.PHONY: distclean
distclean: clean
	-rm -rf $(LIB) .update_header

-include *.d
