include ../config.make

HEADERS = kmrnext.hpp
LIB     = libkmrnext.a
OBJS    = basic_data.o kmrnext.o
TESTS   = test00-basic test01-mng-ds

ifeq ($(TARGET), KMR)
    DEPLIB += $(KMRDIR)/lib/libkmr.a
endif

.Phony: all
all: $(LIB)

.Phony: tests
tests: $(TESTS)

libkmrnext.a: .update_header $(HEADERS) $(OBJS)
	ar csr $@ $(OBJS)

.cpp.o:
ifeq ($(TARGET), KMR)
	$(CXX) $(DEFS) $(CXXFLAGS) -I$(KMRDIR)/include -c $<
else
	$(CXX) $(DEFS) $(CXXFLAGS) -c $<
endif

test00-basic: test00-basic.o $(LIB)
	$(CXX) $(LDFLAGS) -o $@ $^ $(DEPLIB)

test01-mng-ds: test01-mng-ds.o $(LIB)
	$(CXX) $(LDFLAGS) -o $@ $^ $(DEPLIB)

.update_header:
	cp kmrnext.hpp _kmrnext.hpp
	sed -e "s/#define BACKEND_[A-Z]* 1/#define BACKEND_$(TARGET) 1/" \
	< _kmrnext.hpp > kmrnext.hpp
	touch .update_header
	rm -f _kmrnext.hpp

.PHONY: clean
clean:
	-rm -rf $(OBJS) $(TESTS) *.o *.d

.PHONY: distclean
distclean: clean
	-rm -rf $(LIB) .update_header

-include *.d

