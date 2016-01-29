include ../config.make

LIB    = ../src/libkmrnext.a
PROGS  = nicam-letkf-pattern nicam-letkf-pattern-kmr \
	 viz-pattern viz-pattern-kmr \
	 remd-pattern remd-pattern-kmr
DEPLIB = $(KMRDIR)/lib/libkmr.a

.Phony: all
all: $(PROGS)

nicam-letkf-pattern: nicam-letkf-pattern.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIB) $(DEPLIB)

viz-pattern: viz-pattern.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIB) $(DEPLIB)

remd-pattern: remd-pattern.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIB) $(DEPLIB)

nicam-letkf-pattern-kmr: nicam-letkf-pattern-kmr.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIB) $(DEPLIB)

viz-pattern-kmr: viz-pattern-kmr.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIB) $(DEPLIB)

remd-pattern-kmr: remd-pattern-kmr.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIB) $(DEPLIB)

.cpp.o:
ifneq ($(TARGET), KMR)
	$(error KMR backend should be used)
endif
	$(CXX) $(DEFS) $(CXXFLAGS) -I../src -I$(KMRDIR)/include -c $<

.PHONY: clean
clean:
	-rm -rf $(PROGS) *.o *.d

.PHONY: distclean
distclean: clean
