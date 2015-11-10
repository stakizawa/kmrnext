include ../config.make

HEADERS = kmrnext.hpp kmrnext_impl.hpp
LIB     = libkmrnext.a
OBJS    = kmrnext.o basic_data.o data_store.o
TESTS   = test00-basic test01-mng-ds

ifeq ($(TARGET), KMR)
    DEPLIB += $(KMRDIR)/lib/libkmr.a
endif

.Phony: all
all: $(LIB)

.Phony: tests
tests: $(TESTS)

libkmrnext.a: $(HEADERS) $(OBJS)
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

.PHONY: clean
clean:
	-rm -rf $(OBJS) $(TESTS) *.o *.d

.PHONY: distclean
distclean: clean
	-rm -rf $(LIB)

-include *.d

