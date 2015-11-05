include ../config.make

LIB   = ../src/libkmrnext.a
PROGS = nicam-letkf-pattern

ifeq ($(TARGET), KMR)
    DEPLIB += $(KMRDIR)/lib/libkmr.a
endif

.Phony: all
all: $(PROGS)

nicam-letkf-pattern: nicam-letkf-pattern.o
	$(CXX) -o $@ $^ $(LIB) $(DEPLIB)

.cpp.o:
	$(CXX) $(DEFS) $(CXXFLAGS) -I../src -c $<

.PHONY: clean
clean:
	-rm -rf $(PROGS) *.o *.d

.PHONY: distclean
distclean: clean
