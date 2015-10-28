include ../config.make

HEADERS = kmrnext.hpp
OBJS    = kmrnext.o
TESTS   = test00-basic test01-mng-ds

.Phony: all
all: $(TESTS)

test00-basic: $(OBJS) test00-basic.o
	$(CXX) -o $@ $^

test01-mng-ds: $(OBJS) test01-mng-ds.o
	$(CXX) -o $@ $^

.cpp.o:
	$(CXX) $(DEFS) $(CFLAGS) -c $<

.PHONY: clean
clean:
	-rm -rf $(OBJS)
	-rm -rf $(TESTS) *.o *.d

.PHONY: distclean
distclean: clean

-include *.d

