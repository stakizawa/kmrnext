include ../config.make

LIB   = ../src/libkmrnext.a
PROGS = nicam-letkf-pattern in-situ-viz-pattern remd-pattern

ifeq ($(TARGET), KMR)
    DEPLIB += $(KMRDIR)/lib/libkmr.a
endif

.Phony: all
all: $(PROGS)

nicam-letkf-pattern: nicam-letkf-pattern.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIB) $(DEPLIB)

in-situ-viz-pattern: in-situ-viz-pattern.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIB) $(DEPLIB)

remd-pattern: remd-pattern.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIB) $(DEPLIB)

.cpp.o:
ifeq ($(TARGET), KMR)
	$(CXX) $(DEFS) $(CXXFLAGS) -I../src -I$(KMRDIR)/include -c $<
else
	$(CXX) $(DEFS) $(CXXFLAGS) -I../src -c $<
endif

.PHONY: clean
clean:
	-rm -rf $(PROGS) *.o *.d

.PHONY: distclean
distclean: clean
