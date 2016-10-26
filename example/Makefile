include ../config.make

LIB   = ../src/libkmrnext.a
PROGS = nicam-letkf-pattern in-situ-viz-pattern remd-pattern \
	matrix-transpose qcdmd-pattern

ifeq ($(TARGET), KMR)
    DEPLIB += $(KMRDIR)/lib/libkmr.a
    PROGS  += bqcd-pattern
endif

.Phony: all
all: $(PROGS)

nicam-letkf-pattern: nicam-letkf-pattern.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIB) $(DEPLIB)

in-situ-viz-pattern: in-situ-viz-pattern.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIB) $(DEPLIB)

remd-pattern: remd-pattern.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIB) $(DEPLIB)

matrix-transpose: matrix-transpose.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIB) $(DEPLIB)

qcdmd-pattern: qcdmd-pattern.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIB) $(DEPLIB)

bqcd-pattern: bqcd-pattern.o
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
