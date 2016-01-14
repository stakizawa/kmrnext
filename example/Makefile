include ../config.make

LIB   = ../src/libkmrnext.a
PROGS = nicam-letkf-pattern in-situ-viz-pattern

ifeq ($(TARGET), KMR)
    PROGS  += bm-nicam-letkf-pattern bm-in-situ-viz-pattern
    DEPLIB += $(KMRDIR)/lib/libkmr.a
endif

.Phony: all
all: $(PROGS)

nicam-letkf-pattern: nicam-letkf-pattern.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIB) $(DEPLIB)

in-situ-viz-pattern: in-situ-viz-pattern.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIB) $(DEPLIB)


bm-nicam-letkf-pattern: bm-nicam-letkf-pattern.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIB) $(DEPLIB)

bm-in-situ-viz-pattern: bm-in-situ-viz-pattern.o
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
