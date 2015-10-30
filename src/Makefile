include ../config.make

HEADERS = kmrnext.hpp
LIB     = libkmrnext.a
OBJS    = kmrnext.o
TESTS   = test00-basic test01-mng-ds

.Phony: all
all: $(LIB)

.Phony: tests
tests: $(TESTS)

libkmrnext.a: kmrnext.hpp $(OBJS)
	ar csr $@ $(OBJS)

.cpp.o:
	$(CXX) $(DEFS) $(CXXFLAGS) -c $<

test00-basic: test00-basic.o $(LIB)
	$(CXX) -o $@ $^

test01-mng-ds: test01-mng-ds.o $(LIB)
	$(CXX) -o $@ $^

.PHONY: clean
clean:
	-rm -rf $(OBJS) $(TESTS) *.o *.d

.PHONY: distclean
distclean: clean
	-rm -rf $(LIB)

-include *.d

